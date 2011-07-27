'''
Created on Jun 13, 2011
@author: Douglas Crandell
'''
from __future__ import division
from Bio.Seq import Seq
import math
import os
import sqlite3
import sys
import wx
import itertools
import threading
import re
import random
import genome_space as gs
import scipy #@UnresolvedImport
import matplotlib.pyplot as plt #@UnresolvedImport
from mpl_toolkits.mplot3d import Axes3D #@UnresolvedImport
DISTANCE_ITERATIONS=5

class Input():
    def add_sequences(self):
        #Launch dialog for getting FASTA files to add sequences to database
        filters = 'FASTA files (*.fasta; *.fastq; *.fa; *.fna; *.seq; *.txt)|*.fasta;*.fastq;*.fa;*.fna;*.seq;*.txt'
        dialog = wx.FileDialog(None, style = wx.OPEN, wildcard = filters)
        if dialog.ShowModal() == wx.ID_OK:
            self.fasta = dialog.GetPath()
            Database.create_genome_table(Database(), self.fasta)
        else:
            self.fasta = ""
        return self.fasta
    
    def input_file(self):
        #Launch dialog for getting FASTA input file for quick calc
        filters = 'FASTA files (*.fasta; *.fastq; *.fa; *.fna; *.seq; *.txt)|*.fasta;*.fastq;*.fa;*.fna;*.seq;*.txt'
        dialog = wx.FileDialog(None, style = wx.OPEN, wildcard = filters)
        if dialog.ShowModal() == wx.ID_OK:
            self.fasta = dialog.GetPath()
        else:
            self.fasta = ""
        return self.fasta
        
    def read_seqs_from_file(self, filename):
        #Reads sequences from file, returns name of sequence, sequence id, the sequence, and the accession number
        id, name, sequence, accession = [], [], [], []
        seq = ""
        with open(filename) as file:
            for line in file:
                if line.startswith('>'):
                    id.append(line.split(' ')[0].lstrip('>'))
                    if "|" in line:
                        accession.append(line.split('|')[3].rstrip('.'))
                    else:
                        accession.append(line.split(' ')[0].lstrip('>'))
                    name.append(line[line.find(' '):].split(',')[0])
                    if seq != "":
                        sequence.append(seq)
                        seq = ""
                else:
                    seq = seq + line.rstrip('\n')
            sequence.append(seq)
        return name, id, sequence, accession

class Dictionary():
    #Add word to dictionary. If word already in dictionary, then increase word count
    def add_word(self, dict, word, pagenumber):
        dict[word] = dict.get(word, 0) + float(pagenumber)

class Database():
    def generalize_seq(self, sequence):
        #Take nucleotide sequence and generalize bases as purines or pyrimidines.
        #Could possibly allow classification on a broader level by eliminating/reducing strain or species specific differences
        new_seq = ""
        for char in sequence:
            if char == 'A' or char == 'G':
                new_seq += 'R'
            elif char == 'T' or char == 'C':
                new_seq += 'Y'
        return new_seq
    
    def js_range(self, seq1, seq2):
        #Computes Jensen-Shannon Divergence for a range of kmers for k = 3 to k = 8
        F = Frequency()
        dist_list = []
        for k in xrange(3,9):
            dist_list.append(F.jensen_shannon(F.frequency_profile(seq1,k)[1], F.frequency_profile(seq2,k)[1]))
        return dist_list
    
    def frag_seq(self, seq1, seq2):
        frag, i, distances, result = "", 0, [], []
        for char in seq1:
            if i < len(seq2):
                frag += char
                i += 1
            else:
                frag += char
                dists = self.js_range(frag, seq2)
                distances.append(dists)
                i = 0
                frag = ""
        dists = self.js_range(frag, seq2)
        distances.append(dists)
        total, sum = 0, 0
        for j in range(0,6):
            for list in distances:
                total += 1
                sum += list[j]
            result.append(sum/total)
        return result
    
    def vector_range(self, sequence, type):
        #Computes moment vector for a sequence and natural vector for a range of kmers, k=1 to k = 15
        V = gs.Vectors()
        fp_list, v_list = [], []
        if type == True:
            ycoords = V.nt_graph(sequence)
        elif type == False:
            ycoords = V.aa_graph(sequence)
        ycoord_string = ','.join(str(n) for n in ycoords)
        v_list.append(ycoord_string)
        mv = V.moment_vector(ycoords)
        moment_string = ','.join(str(n) for n in mv)
        v_list.append(moment_string)
        if type == False:
            sequence = str(Seq(sequence).translate())
        for k in xrange(1,6):
            fp_list.append(V.frequency_profile(sequence, k))
        for profile in fp_list:
            nv = V.natural_vector(sequence, profile)
            nv_string = ','.join(str(n) for n in nv)
            v_list.append(nv_string)
        return v_list
    
    def distmatrix(self, id_list, distances):
        #Create distance matrix
        matrix = scipy.zeros((len(id_list), len(id_list)))
        num = 0
        for id in id_list:
            n, i, d = 0, 0, [] 
            for dist in distances[id]:
                if str(dist) == '0.00' and i == 0:
                    d.append(dist)
                    i += 1
                elif str(dist) != '0.00':
                    d.append(dist)
                else:
                    pass
            for dist in d:
                matrix[num,n] = dist
                n +=1
            num += 1
        return matrix
    
    def fastmap(self, dist, K):
        #dist is a NxN distance matrix returns coordinates for each N in K dimensions
        return FastMap(dist,True).map(K)
    
    def vlen(self,v):
        return math.sqrt(scipy.vdot(v,v))
        
    def create_genome_table(self, filename):
        #Import sequences, accession #'s, and sequence names into database
        #Currently also import sequence id, probably will be removed in the future
        accession_list, seq_list = [], []
        reads = Input().read_seqs_from_file(filename)
        connection = sqlite3.connect('distance.db')
        cursor = connection.cursor()
        cursor.execute('''CREATE TABLE IF NOT EXISTS genomes (id integer primary key, sequence_name, accession UNIQUE, sequence_id, nt_sequence)''')
        for seq_name, seq_id, nt_seq, accession in itertools.izip(reads[0], reads[1], reads[2], reads[3]):
            t = seq_name, accession, seq_id, nt_seq
            prev = cursor.execute('''SELECT COUNT(accession) FROM genomes''')
            for row in prev:
                prev = row[0]
            cursor.execute('''INSERT OR IGNORE INTO genomes (id, sequence_name, accession, sequence_id, nt_sequence) VALUES (NULL, ?,?,?,?)''',t)
            current = cursor.execute('''SELECT COUNT(accession) FROM genomes''')
            for row in current:
                current = row[0]
            if prev != current:
                accession_list.append(accession)
                seq_list.append(nt_seq)
        connection.commit()
        self.create_distance_tables(accession_list,seq_list)
        self.create_js_table(accession_list,seq_list)
        self.create_vector_tables(accession_list, seq_list)
        
    def create_js_table(self, ids, sequences):
        #Import Jensen-Shannon divergence distances into database
        old_ids, old_sequences = [], []
        old_ids.extend(ids)
        connection = sqlite3.connect('distance.db')
        cursor = connection.cursor()
        cursor.execute('''CREATE TABLE IF NOT EXISTS nt_js (id integer primary key, from_accession, to_accession, three, four, five, six, seven, eight)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS aa_js (id integer primary key, from_accession, to_accession, three, four, five, six, seven, eight)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS gen_js (id integer primary key, from_accession, to_accession, three, four, five, six, seven, eight)''')
        cursor.execute('''SELECT accession, nt_sequence FROM genomes''')
        
        for row in cursor:
            old_ids.append(row[0])
            old_sequences.append(row[1])
        nt_org_dict, aa_org_dict, gen_org_dict = {}, {}, {}
        num = 0
        total = len(sequences)*len(old_sequences)
        dialog = wx.ProgressDialog("Jensen-Shannon Calculations", "Time Remaining", total, style=wx.PD_ELAPSED_TIME|wx.PD_REMAINING_TIME)
        for id1, seq1 in itertools.izip(ids, sequences):
            nt_dict, aa_dict, gen_dict = {}, {}, {}
            dialog.Update(num)
            aa_seq1 = str(Seq(seq1).translate())
            gen_seq1 = self.generalize_seq(seq1)
            for id2, seq2 in itertools.izip(old_ids, old_sequences):
                num += 1
                dialog.Update(num)
                aa_seq2 = str(Seq(seq2).translate())
                gen_seq2 = self.generalize_seq(seq2)
                if id2 not in nt_org_dict:
                    if id1 == id2:
                        nt_dict[id2] = [0,0,0,0,0,0]
                        aa_dict[id2] = [0,0,0,0,0,0]
                        gen_dict[id2] = [0,0,0,0,0,0]
                    else:
                        if len(seq1) > 3*len(seq2):
                            nt_dict[id2] = self.frag_seq(seq1, seq2)
                            aa_dict[id2] = self.frag_seq(aa_seq1, aa_seq2)
                            gen_dict[id2] = self.frag_seq(gen_seq1, gen_seq2)
                        elif len(seq2) > 3*len(seq1):
                            nt_dict[id2] = self.frag_seq(seq2, seq1)
                            aa_dict[id2] = self.frag_seq(aa_seq2, aa_seq1)
                            gen_dict[id2] = self.frag_seq(gen_seq2, gen_seq1)
                        else:
                            nt_dict[id2] = self.js_range(seq1, seq2)
                            aa_dict[id2] = self.js_range(aa_seq1, aa_seq2)
                            gen_dict[id2] = self.js_range(gen_seq1, gen_seq2)
                else:
                    pass
            nt_org_dict[id1], aa_org_dict[id1], gen_org_dict[id1] = nt_dict, aa_dict, gen_dict
                
        for d1, d2, d3 in zip(nt_org_dict.iteritems(),aa_org_dict.iteritems(),gen_org_dict.iteritems()):
            org = d1[0]
            for dest, val in d1[1].iteritems():
                a = org, dest, val[0], val[1], val[2], val[3], val[4], val[5]
                cursor.execute('''INSERT INTO nt_js (id, from_accession, to_accession, three, four, five, six, seven, eight) VALUES (NULL,?,?,?,?,?,?,?,?)''',a)
            for dest, val in d2[1].iteritems():
                a = org, dest, val[0], val[1], val[2], val[3], val[4], val[5]
                cursor.execute('''INSERT INTO aa_js (id, from_accession, to_accession, three, four, five, six, seven, eight) VALUES (NULL,?,?,?,?,?,?,?,?)''',a)
            for dest, val in d3[1].iteritems():
                a = org, dest, val[0], val[1], val[2], val[3], val[4], val[5]
                cursor.execute('''INSERT INTO gen_js (id, from_accession, to_accession, three, four, five, six, seven, eight) VALUES (NULL,?,?,?,?,?,?,?,?)''',a)
        connection.commit()
        dialog.Destroy()
        
    def create_vector_tables(self, ids, sequences):
        #Import vectors and graphs into database
        old_ids, num =  [], 0
        connection = sqlite3.connect('distance.db')
        cursor = connection.cursor()
        cursor.execute('''CREATE TABLE IF NOT EXISTS nt_vectors (id integer primary key, accession, ycoords, moment, one, two, three, four, five)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS aa_vectors (id integer primary key, accession, ycoords, moment, one, two, three, four, five)''')
        cursor.execute('''SELECT accession FROM nt_vectors''')
        for row in cursor:
            old_ids.append(row[0])
        dialog = wx.ProgressDialog("Vector Calculations", "Time Remaining", len(ids), style=wx.PD_ELAPSED_TIME|wx.PD_REMAINING_TIME)
        for id, sequence in itertools.izip(ids,sequences):
            num += 1
            dialog.Update(num)
            if id not in old_ids:
                nt_vectors = self.vector_range(sequence, True)
                aa_vectors = self.vector_range(sequence, False)
                nt = id, nt_vectors[0], nt_vectors[1], nt_vectors[2], nt_vectors[3], nt_vectors[4], nt_vectors[5], nt_vectors[6]
                aa = id, aa_vectors[0], aa_vectors[1], aa_vectors[2], aa_vectors[3], aa_vectors[4], aa_vectors[5], aa_vectors[6]
                cursor.execute('''INSERT INTO nt_vectors (id, accession, ycoords, moment, one, two, three, four, five) VALUES (NULL,?,?,?,?,?,?,?,?)''',nt)
                cursor.execute('''INSERT INTO aa_vectors (id, accession, ycoords, moment, one, two, three, four, five) VALUES (NULL,?,?,?,?,?,?,?,?)''',aa)
        connection.commit()
        dialog.Destroy()
        dialog = wx.MessageDialog(None, "%s sequences added to database" %num, "Import Complete",wx.OK)
        dialog.ShowModal()
        dialog.Destroy()

    def create_distance_tables(self, ids, sequences):
        #Calculate distances for frequency distribution calculations and store in database
        F = Frequency()
        nt_org_dict, aa_org_dict, gen_org_dict = {}, {}, {}
        old_ids, old_sequences = [], []
        old_ids.extend(ids)
        num = 0

        connection = sqlite3.connect('distance.db')
        cursor = connection.cursor()
        cursor.execute('''CREATE TABLE IF NOT EXISTS nt_odds_ratio (id integer primary key, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS nt_log_odds (id integer primary key, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS nt_odds_diff (id integer primary key, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS nt_poisson (id integer primary key, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS aa_odds_ratio (id integer primary key, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS aa_log_odds (id integer primary key, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS aa_odds_diff (id integer primary key,from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS aa_poisson (id integer primary key, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS gen_odds_ratio (id integer primary key, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS gen_log_odds (id integer primary key, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS gen_odds_diff (id integer primary key,from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS gen_poisson (id integer primary key, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist)''')
        cursor.execute('''SELECT accession, nt_sequence FROM genomes''')
        
        for row in cursor:
            old_ids.append(row[0])
            old_sequences.append(row[1])
        total = len(ids)*len(old_sequences)
        
        dialog = wx.ProgressDialog("K-mer Calculations", "Time Remaining", total, style=wx.PD_ELAPSED_TIME|wx.PD_REMAINING_TIME)
        for id1, seq1, in itertools.izip(ids, sequences):
            nt_dict, aa_dict, gen_dict = {}, {}, {}
            aa_seq1 = str(Seq(seq1).translate())
            gen_seq1 = self.generalize_seq(seq1)
            dialog.Update(num)
            for id2, seq2 in itertools.izip(old_ids, old_sequences):
                num += 1
                dialog.Update(num)
                aa_seq2 = str(Seq(seq2).translate())
                gen_seq2 = self.generalize_seq(seq2)
                if id2 not in nt_org_dict:
                    if id1 == id2:
                        nt_dict[id2] = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
                        aa_dict[id2] = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
                        gen_dict[id2] = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
                    else:
                        nt_dict[id2] = [F.odds(2,seq1,seq2),F.odds(3,seq1,seq2),F.odds(4,seq1,seq2),F.odds(5,seq1,seq2)]
                        aa_dict[id2] = [F.odds(2,aa_seq1,aa_seq2),F.odds(3,aa_seq1,aa_seq2),F.odds(4,aa_seq1,aa_seq2),F.odds(5,aa_seq1,aa_seq2)]
                        gen_dict[id2] = [F.odds(2,gen_seq1,gen_seq2),F.odds(3,gen_seq1,gen_seq2),F.odds(4,gen_seq1,gen_seq2),F.odds(5,gen_seq1,gen_seq2)]
                else:
                    pass
            nt_org_dict[id1], aa_org_dict[id1], gen_org_dict[id1] = nt_dict, aa_dict, gen_dict
                    
        for d1, d2, d3 in zip(nt_org_dict.iteritems(),aa_org_dict.iteritems(),gen_org_dict.iteritems()):
            org = d1[0]
            for dest, val in d1[1].iteritems():
                a = org, dest, val[0][0], val[1][0], val[2][0], val[3][0]
                b = org, dest, val[0][1], val[1][1], val[2][1], val[3][1]
                c = org, dest, val[0][2], val[1][2], val[2][2], val[3][2]
                d = org, dest, val[0][3], val[1][3], val[2][3], val[3][3]
                cursor.execute('''INSERT INTO nt_odds_ratio (id, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist) VALUES (NULL, ?,?,?,?,?,?)''', a)
                cursor.execute('''INSERT INTO nt_log_odds (id, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist) VALUES (NULL, ?,?,?,?,?,?)''', b)
                cursor.execute('''INSERT INTO nt_odds_diff (id, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist) VALUES (NULL, ?,?,?,?,?,?)''', c)
                cursor.execute('''INSERT INTO nt_poisson (id, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist) VALUES (NULL, ?,?,?,?,?,?)''', d)
            for dest, val in d2[1].iteritems():
                a = org, dest, val[0][0], val[1][0], val[2][0], val[3][0]
                b = org, dest, val[0][1], val[1][1], val[2][1], val[3][1]
                c = org, dest, val[0][2], val[1][2], val[2][2], val[3][2]
                d = org, dest, val[0][3], val[1][3], val[2][3], val[3][3]
                cursor.execute('''INSERT INTO aa_odds_ratio (id, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist) VALUES (NULL, ?, ?, ?, ?, ?, ?)''',a)
                cursor.execute('''INSERT INTO aa_log_odds (id, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist) VALUES (NULL, ?, ?, ?, ?, ?, ?)''', b)
                cursor.execute('''INSERT INTO aa_odds_diff (id, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist) VALUES (NULL, ?, ?, ?, ?, ?, ?)''',c)
                cursor.execute('''INSERT INTO aa_poisson (id, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist) VALUES (NULL, ?, ?, ?, ?, ?, ?)''',d)
            for dest, val in d3[1].iteritems():
                a = org, dest, val[0][0], val[1][0], val[2][0], val[3][0]
                b = org, dest, val[0][1], val[1][1], val[2][1], val[3][1]
                c = org, dest, val[0][2], val[1][2], val[2][2], val[3][2]
                d = org, dest, val[0][3], val[1][3], val[2][3], val[3][3]
                cursor.execute('''INSERT INTO gen_odds_ratio (id, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist) VALUES (NULL, ?, ?, ?, ?, ?, ?)''',a)
                cursor.execute('''INSERT INTO gen_log_odds (id, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist) VALUES (NULL, ?, ?, ?, ?, ?, ?)''', b)
                cursor.execute('''INSERT INTO gen_odds_diff (id, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist) VALUES (NULL, ?, ?, ?, ?, ?, ?)''',c)
                cursor.execute('''INSERT INTO gen_poisson (id, from_accession, to_accession, di_dist, tri_dist, tet_dist, pent_dist) VALUES (NULL, ?, ?, ?, ?, ?, ?)''',d)
        connection.commit()
        dialog.Destroy()
    
    def quick_calc(self, k, organisms_list, type, calculation):
        #Method to draw tree from FASTA file for ONE selected tree and calculation type quickly without importing sequences into database and computing
        #other calculations
        distances, done = {}, []
        if calculation != 'moment' and calculation != 'natural':
            F = Frequency()
            for accession1, seq1 in itertools.izip(organisms_list[3], organisms_list[2]):
                if type == 'aa':
                    seq1 = str(Seq(seq1).translate())
                for accession2, seq2 in itertools.izip(organisms_list[3], organisms_list[2]):
                    if type == 'aa':
                        seq2 = str(Seq(seq2).translate())
                    if accession1 == accession2:
                        if accession1 not in distances:
                            distances[accession1] = [0.00]
                        else:
                            tmp = distances[accession1]
                            tmp.append(0.00)
                            distances[accession1] = tmp
                    elif accession2 not in done:
                        if calculation == 'oddsratio':
                            tmp = distances[accession1]
                            tmp.append(F.odds(k, seq1, seq2)[0])
                            distances[accession1] = tmp
                            if accession2 not in distances:
                                distances[accession2] = [F.odds(k, seq1, seq2)[0]]
                            else:
                                tmp = distances[accession2]
                                tmp.append(F.odds(k, seq1, seq2)[0])
                                distances[accession2] = tmp
                        elif calculation == 'logodds':
                            tmp = distances[accession1]
                            tmp.append(F.odds(k, seq1, seq2)[1])
                            distances[accession1] = tmp
                            if accession2 not in distances:
                                distances[accession2] = [F.odds(k, seq1, seq2)[1]]
                            else:
                                tmp = distances[accession2]
                                tmp.append(F.odds(k, seq1, seq2)[1])
                                distances[accession2] = tmp
                        elif calculation == 'oddsdiff':
                            tmp = distances[accession1]
                            tmp.append(F.odds(k, seq1, seq2)[2])
                            distances[accession1] = tmp
                            if accession2 not in distances:
                                distances[accession2] = [F.odds(k, seq1, seq2)[2]]
                            else:
                                tmp = distances[accession2]
                                tmp.append(F.odds(k, seq1, seq2)[2])
                                distances[accession2] = tmp
                        elif calculation == 'poisson':
                            tmp = distances[accession1]
                            tmp.append(F.odds(k, seq1, seq2)[3])
                            distances[accession1] = tmp
                            if accession2 not in distances:
                                distances[accession2] = [F.odds(k, seq1, seq2)[3]]
                            else:
                                tmp = distances[accession2]
                                tmp.append(F.odds(k, seq1, seq2)[3])
                                distances[accession2] = tmp
                        elif calculation == 'jensen_shannon':
                            tmp = distances[accession1]
                            tmp.append(F.jensen_shannon(F.frequency_profile(seq1,k)[0], F.frequency_profile(seq2,k)[0]))
                            distances[accession1] = tmp
                            if accession2 not in distances:
                                distances[accession2] = [F.jensen_shannon(F.frequency_profile(seq1,k)[0], F.frequency_profile(seq2,k)[0])]
                            else:
                                tmp = distances[accession2]
                                tmp.append(F.jensen_shannon(F.frequency_profile(seq1,k)[0], F.frequency_profile(seq2,k)[0]))
                                distances[accession2] = tmp
                done.append(accession1)
        else:
            V = gs.Vectors()
            v_dict, distances = {}, {}
            for acc, seq in itertools.izip(organisms_list[3], organisms_list[2]):
                if calculation == 'moment':
                    if type == 'nt':
                        mv = V.moment_vector(V.nt_graph(seq))
                    else:
                        mv = V.moment_vector(V.aa_graph(seq))
                    moment_string = ','.join(str(n) for n in mv)
                    v_dict[acc] = moment_string
                elif calculation == 'natural':
                    if type == 'nt':
                        nv = V.natural_vector(seq, V.frequency_profile(seq, k))
                    else:
                        seq = str(Seq(seq).translate())
                        nv = V.natural_vector(seq, V.frequency_profile(seq, k))
                    natural_string = ','.join(str(n) for n in nv)
                    v_dict[acc] = natural_string
            for v1 in organisms_list[3]:
                for v2 in organisms_list[3]:
                    if v1 == v2:
                        if v1 not in distances:
                            distances[v1] = [0.00]
                        else:
                            tmp = distances[v1]
                            tmp.append(0.00)
                            distances[v1] = tmp
                    elif v2 not in done:
                        tmp = distances[v1]
                        tmp.append(V.vector_distance(eval(v_dict[v1]), eval(v_dict[v2])))
                        distances[v1] = tmp
                        if v2 not in distances:
                            distances[v2] = [V.vector_distance(eval(v_dict[v1]), eval(v_dict[v2]))]
                        else:
                            tmp = distances[v2]
                            tmp.append(V.vector_distance(eval(v_dict[v1]), eval(v_dict[v2])))
                            distances[v2] = tmp       
                done.append(v1)
                    
        file = open('infile', 'w')
        file.write(str(len(organisms_list[3])) + '\n')
        for org in organisms_list[3]:
            string = ""
            for d in distances[org]:
                d = "%.2f" %d
                string += str(d) + ' '
            string.strip(' ')
            org = org.strip(' ')
            file.write(org[:9] + ' '+ string + '\n')
        file.close()
            
    def retrieve_distances(self, k, organisms_list, type, calculation, bool):
        #Retrieve distances from database, create distance matrix, and draw tree
        connection = sqlite3.connect('distance.db')
        cursor = connection.cursor()
        if calculation != 'moment' and calculation != 'natural':
            if type == 'nucleotide':
                if calculation == 'oddsratio':
                    if k == 2:
                        cursor.execute('''SELECT from_accession,to_accession,di_dist FROM nt_odds_ratio''')
                    elif k ==3:
                        cursor.execute('''SELECT from_accession,to_accession,tri_dist FROM nt_odds_ratio''')
                    elif k==4:
                        cursor.execute('''SELECT from_accession,to_accession,tet_dist FROM nt_odds_ratio''')
                    elif k==5:
                        cursor.execute('''SELECT from_accession,to_accession,pent_dist FROM nt_odds_ratio''')
                elif calculation == 'logodds':
                    if k == 2:
                        cursor.execute('''SELECT from_accession,to_accession,di_dist FROM nt_log_odds''')
                    elif k ==3:
                        cursor.execute('''SELECT from_accession,to_accession,tri_dist FROM nt_log_odds''')
                    elif k==4:
                        cursor.execute('''SELECT from_accession,to_accession,tet_dist FROM nt_log_odds''')
                    elif k==5:
                        cursor.execute('''SELECT from_accession, to_accession, pent_dist FROM nt_log_odds''')
                elif calculation == 'oddsdiff':
                    if k == 2:
                        cursor.execute('''SELECT from_accession,to_accession,di_dist FROM nt_odds_diff''')
                    elif k ==3:
                        cursor.execute('''SELECT from_accession,to_accession,tri_dist FROM nt_odds_diff''')
                    elif k==4:
                        cursor.execute('''SELECT from_accession,to_accession,tet_dist FROM nt_odds_diff''')
                    elif k==5:
                        cursor.execute('''SELECT from_accession, to_accession, pent_dist FROM nt_odds_diff''')
                elif calculation == 'poisson':
                    if k == 2:
                        cursor.execute('''SELECT from_accession,to_accession,di_dist FROM nt_poisson''')
                    elif k ==3:
                        cursor.execute('''SELECT from_accession,to_accession,tri_dist FROM nt_poisson''')
                    elif k==4:
                        cursor.execute('''SELECT from_accession,to_accession,tet_dist FROM nt_poisson''')
                    elif k==5:
                        cursor.execute('''SELECT from_accession, to_accession, pent_dist FROM nt_poisson''')
                elif calculation == 'jensen_shannon':
                    if k == 3:
                        cursor.execute('''SELECT from_accession, to_accession, three FROM nt_js''')
                    elif k ==4:
                        cursor.execute('''SELECT from_accession, to_accession, four FROM nt_js''')
                    elif k == 5:
                        cursor.execute('''SELECT from_accession, to_accession, five FROM nt_js''')
                    elif k == 6:
                        cursor.execute('''SELECT from_accession, to_accession, six FROM nt_js''')
                    elif k == 7:
                        cursor.execute('''SELECT from_accession, to_accession, seven FROM nt_js''')
                    elif k == 8:
                        cursor.execute('''SELECT from_accession, to_accession, eight FROM nt_js''')
            elif type == 'aminoacid':
                if calculation == 'oddsratio':
                    if k == 2:
                        cursor.execute('''SELECT from_accession,to_accession,di_dist FROM aa_odds_ratio''')
                    elif k ==3:
                        cursor.execute('''SELECT from_accession,to_accession,tri_dist FROM aa_odds_ratio''')
                    elif k==4:
                        cursor.execute('''SELECT from_accession,to_accession,tet_dist FROM aa_odds_ratio''')
                    elif k==5:
                        cursor.execute('''SELECT from_accession,to_accession,pent_dist FROM aa_odds_ratio''')
                elif calculation == 'logodds':
                    if k == 2:
                        cursor.execute('''SELECT from_accession,to_accession,di_dist FROM aa_log_odds''')
                    elif k ==3:
                        cursor.execute('''SELECT from_accession,to_accession,tri_dist FROM aa_log_odds''')
                    elif k==4:
                        cursor.execute('''SELECT from_accession,to_accession,tet_dist FROM aa_log_odds''')
                    elif k==5:
                        cursor.execute('''SELECT from_accession, to_accession, pent_dist FROM aa_log_odds''')
                elif calculation == 'oddsdiff':
                    if k == 2:
                        cursor.execute('''SELECT from_accession,to_accession,di_dist FROM aa_odds_diff''')
                    elif k ==3:
                        cursor.execute('''SELECT from_accession,to_accession,tri_dist FROM aa_odds_diff''')
                    elif k==4:
                        cursor.execute('''SELECT from_accession,to_accession,tet_dist FROM aa_odds_diff''')
                    elif k==5:
                        cursor.execute('''SELECT from_accession, to_accession, pent_dist FROM aa_odds_diff''')
                elif calculation == 'poisson':
                    if k == 2:
                        cursor.execute('''SELECT from_accession,to_accession,di_dist FROM aa_poisson''')
                    elif k ==3:
                        cursor.execute('''SELECT from_accession,to_accession,tri_dist FROM aa_poisson''')
                    elif k==4:
                        cursor.execute('''SELECT from_accession,to_accession,tet_dist FROM aa_poisson''')
                    elif k==5:
                        cursor.execute('''SELECT from_accession, to_accession, pent_dist FROM nt_poisson''')
                elif calculation == 'jensen_shannon':
                    if k == 3:
                        cursor.execute('''SELECT from_accession, to_accession, three FROM aa_js''')
                    elif k ==4:
                        cursor.execute('''SELECT from_accession, to_accession, four FROM aa_js''')
                    elif k == 5:
                        cursor.execute('''SELECT from_accession, to_accession, five FROM aa_js''')
                    elif k == 6:
                        cursor.execute('''SELECT from_accession, to_accession, six FROM aa_js''')
                    elif k == 7:
                        cursor.execute('''SELECT from_accession, to_accession, seven FROM aa_js''')
                    elif k == 8:
                        cursor.execute('''SELECT from_accession, to_accession, eight FROM aa_js''')
            elif type == 'generalized':
                if calculation == 'oddsratio':
                    if k == 2:
                        cursor.execute('''SELECT from_accession,to_accession,di_dist FROM gen_odds_ratio''')
                    elif k ==3:
                        cursor.execute('''SELECT from_accession,to_accession,tri_dist FROM gen_odds_ratio''')
                    elif k==4:
                        cursor.execute('''SELECT from_accession,to_accession,tet_dist FROM gen_odds_ratio''')
                    elif k==5:
                        cursor.execute('''SELECT from_accession,to_accession,pent_dist FROM gen_odds_ratio''')
                elif calculation == 'logodds':
                    if k == 2:
                        cursor.execute('''SELECT from_accession,to_accession,di_dist FROM gen_log_odds''')
                    elif k ==3:
                        cursor.execute('''SELECT from_accession,to_accession,tri_dist FROM gen_log_odds''')
                    elif k==4:
                        cursor.execute('''SELECT from_accession,to_accession,tet_dist FROM gen_log_odds''')
                    elif k==5:
                        cursor.execute('''SELECT from_accession, to_accession, pent_dist FROM gen_log_odds''')
                elif calculation == 'oddsdiff':
                    if k == 2:
                        cursor.execute('''SELECT from_accession,to_accession,di_dist FROM gen_odds_diff''')
                    elif k ==3:
                        cursor.execute('''SELECT from_accession,to_accession,tri_dist FROM gen_odds_diff''')
                    elif k==4:
                        cursor.execute('''SELECT from_accession,to_accession,tet_dist FROM gen_odds_diff''')
                    elif k==5:
                        cursor.execute('''SELECT from_accession, to_accession, pent_dist FROM gen_odds_diff''')
                elif calculation == 'poisson':
                    if k == 2:
                        cursor.execute('''SELECT from_accession,to_accession,di_dist FROM gen_poisson''')
                    elif k ==3:
                        cursor.execute('''SELECT from_accession,to_accession,tri_dist FROM gen_poisson''')
                    elif k==4:
                        cursor.execute('''SELECT from_accession,to_accession,tet_dist FROM gen_poisson''')
                    elif k==5:
                        cursor.execute('''SELECT from_accession, to_accession, pent_dist FROM gen_poisson''')
                elif calculation == 'jensen_shannon':
                    if k == 3:
                        cursor.execute('''SELECT from_accession, to_accession, three FROM gen_js''')
                    elif k ==4:
                        cursor.execute('''SELECT from_accession, to_accession, four FROM gen_js''')
                    elif k == 5:
                        cursor.execute('''SELECT from_accession, to_accession, five FROM gen_js''')
                    elif k == 6:
                        cursor.execute('''SELECT from_accession, to_accession, six FROM gen_js''')
                    elif k == 7:
                        cursor.execute('''SELECT from_accession, to_accession, seven FROM gen_js''')
                    elif k == 8:
                        cursor.execute('''SELECT from_accession, to_accession, eight FROM gen_js''')   
            id_list, start, dest, dist, distances = [], [], [], [], {}
            for row in cursor:
                start.append(row[0])
                dest.append(row[1])
                dist.append(row[2])
                start.append(row[1])
                dest.append(row[0])
                dist.append(row[2])
            comp = sorted(zip(start, dest, dist))
            for row in comp:
                if row[0] in organisms_list and row[1] in organisms_list:
                    if row[0] not in id_list:
                        id_list.append(row[0])
                    if row[0] in distances:
                        tmp = distances[row[0]]
                        tmp.append("%.2f" %row[2])
                        distances[row[0]] = tmp
                    else:
                        distances[row[0]] = ["%.2f" %row[2]]
            if bool == True:
                with open('infile', 'w') as file:
                    file.write(str(len(distances)) +'\n')
                    string = ""
                    num = 1
                    for id in id_list:
                        if len(id) < 9:
                            string = id +' '*(10-len(id))
                        else:
                            string  = id[:9] + ' '
                        number = 0
                        for dist in distances[id]:
                            number += 1
                            if number != num:
                                string = string + str(dist) + ' '
                        string = string + '\n'
                        file.write(string)
                        num += 1
            else:
                return id_list, distances
        else:
            if type == 'nucleotide':
                if calculation == 'moment':
                    cursor.execute('''SELECT accession, moment FROM nt_vectors''')
                elif calculation == 'natural':
                    if k == 1:
                        cursor.execute('''SELECT accession, one FROM nt_vectors''')
                    elif k == 2:
                        cursor.execute('''SELECT accession, two FROM nt_vectors''')
                    elif k == 3:
                        cursor.execute('''SELECT accession, three FROM nt_vectors''')
                    elif k == 4:
                        cursor.execute('''SELECT accession, four FROM nt_vectors''')
                    elif k == 5:
                        cursor.execute('''SELECT accession, five FROM nt_vectors''')
            elif type == 'aminoacid':
                if calculation == 'moment':
                    cursor.execute('''SELECT accession, moment FROM aa_vectors''')
                elif calculation == 'natural':
                    if k == 1:
                        cursor.execute('''SELECT accession, one FROM aa_vectors''')
                    elif k == 2:
                        cursor.execute('''SELECT accession, two FROM aa_vectors''')
                    elif k == 3:
                        cursor.execute('''SELECT accession, three FROM aa_vectors''')
                    elif k == 4:
                        cursor.execute('''SELECT accession, four FROM aa_vectors''')
                    elif k == 5:
                        cursor.execute('''SELECT accession, five FROM aa_vectors''')
            vector_dict, distances = {}, {}
            id_list, start, dest, dist, done = [], [], [], [], []
            V = gs.Vectors()
            for row in cursor:
                if row[0] in organisms_list:
                    vector_dict[row[0]] = [n for n in row[1].split(',')]
            for v1 in vector_dict.iterkeys():
                for v2 in vector_dict.iterkeys():
                    if v1 != v2 and v2 not in done:
                        start.append(v1)
                        dest.append(v2)
                        distance = V.vector_distance(vector_dict[v1], vector_dict[v2])
                        dist.append(distance)
                        start.append(v2)
                        dest.append(v1)
                        dist.append(distance)
                done.append(v1)
            comp = sorted(zip(start, dest, dist))
            for row in comp:
                if row[0] not in id_list:
                    id_list.append(row[0])
                if row[0] in distances:
                    tmp = distances[row[0]]
                    tmp.append("%.2f" %row[2])
                    distances[row[0]] = tmp
                else:
                    distances[row[0]] = ["%.2f" %row[2]]
            if bool == True:
                with open('infile', 'w') as file:
                    file.write(str(len(distances)) +'\n')
                    string = ""
                    num = 0
                    for id in id_list:
                        tmp = distances[id]
                        tmp.insert(num, "%.2f" %0.00)
                        if len(id) < 9:
                            string = id +' '*(10-len(id))
                        else:
                            string  = id[:9] + ' '
                        number = 0
                        n = 0
                        for dist in distances[id]:
                            string = string + str(dist) + ' '
                            n += 1
                        string = string + '\n'
                        file.write(string)
                        num += 1
            else:
                return id_list, distances
    
    def graph(self, id_list, distances):
        #Distances is a list of distances from one of the frequency calculations
        #ID_list is a list of organisms selected from the list control
        #Creates distance matrix for the selected organisms and then creates a 2D and 3D plot of the organisms in the space
        name_list = []
        connection = sqlite3.connect('distance.db')
        cursor = connection.cursor()
        cursor.execute('''SELECT sequence_name, accession FROM genomes''')
        for row in cursor:
            if row[1] in id_list:
                name_list.append(row[0])
        dm = self.distmatrix(id_list, distances) #Create distance matrix
        p = self.fastmap(dm, 3)
        plt.scatter([x[0] for x in p], [x[1] for x in p], c='bgrcmyk')
        af = gs.AnnoteFinder([x[0] for x in p], [x[1] for x in p],name_list)
        plt.connect('button_press_event', af) #Connect plot to button press event for annotating points 2D plot on click
        font = {'fontsize':12}
        plt.title('2D Distance Matrix Representation', font)
        plt.xlabel('X')
        plt.ylabel('Y')
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter3D([x[0] for x in p], [x[1] for x in p], map(self.vlen,p), c='red')
        num = 0
        for x,y,z in zip([x[0] for x in p], [x[1] for x in p], map(self.vlen,p)): #Annotate 3D plot
            ax.text3D(x,y,z, name_list[num], fontsize = 8)
            num += 1
        ax.set_title('3D Distance Matrix Representation', font)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()
                        
class Frequency():
    def frequency_profile(self, sequence, k):
        #Compute a frequence distribution profile for kmers in a sequence
        k = int(k)
        length = len(sequence)
        total = length - k + 1
        count = 0
        string = ""
        profile, distribution = {}, {}
        for letter in sequence:
            count += 1
            if count < k:
                string = string + letter
            elif count == k:
                string = string + letter
                Dictionary().add_word(profile, string, 1)
            else:
                string = string[1:] + letter
                Dictionary().add_word(profile, string, 1)
        for word in profile:
            distribution[word] = profile.get(word)/total
        return profile, distribution

    def kldiv(self, p, q):
        #Compute Kullback-Leibler Divergence
        if len(p) == 0:
            return 1e33
        if len(q) == 0:
            return 1e33

        psum = 0 + sum(p.values())
        qsum = 0 + sum(q.values())
        
        vocabdiff = set(p.keys()).difference(set(q.keys()))
        lenvocabdiff = len(vocabdiff)
        epsilon = min(min(p.values())/psum, min(q.values())/qsum) * 0.001
        gamma = 1 - lenvocabdiff * epsilon

        sc = sum([v/psum for v in p.itervalues()])
        st = sum([v/qsum for v in q.itervalues()])
        if sc < 9e-6:
            sys.exit(2)
        if st < 9e-6:
            sys.exit(2)

        div = 0
        for t, v in p.iteritems():
            ptp = v /psum
            ptq = epsilon
            if t in q:
                ptq = gamma * (q[t]/qsum)
            ckl = (ptp - ptq) * math.log((ptp/ptq),2)
            div += ckl
        return div
                
    def jensen_shannon(self, p, q):
        #Compute Jensen Shannon divergence using the Kullback-Leibler divergence to measure the similarity
        #between two probability distributions
        return math.sqrt(0.5*self.kldiv(p,q) + 0.5*self.kldiv(q,p))

    def odds_ratio(self, multi_freq, frequency):
        #Compute bias vector as odds ratio
        nt_bias = {}
        for combo in multi_freq.iterkeys():
            exp_prob = 1
            for letter in combo:
                exp_prob = exp_prob * frequency[letter]
            if exp_prob != 0:
                odds = multi_freq[combo]/exp_prob
            else:
                odds = 0
            Dictionary().add_word(nt_bias, combo, odds)
        return nt_bias

    def log_odds_ratio(self, multi_freq, frequency):
        #Compute bias vector as logarithm of odds ratio
        nt_bias = {}
        for combo in multi_freq.iterkeys():
            exp_prob = 1
            for letter in combo:
                exp_prob = exp_prob * frequency[letter]
            if multi_freq[combo] != 0 and exp_prob != 0:
                odds = math.log10(multi_freq[combo]/exp_prob)
            else:
                odds = 0
            Dictionary().add_word(nt_bias, combo, odds)
        return nt_bias

    def odds_difference(self, multi_freq, frequency):
        #Computes nucleotide bias vector as the odds difference.  Similar to Poisson distribution without the normalization factor
        nt_bias = {}
        for combo in multi_freq.iterkeys():
            exp_prob = 1
            for letter in combo:
                exp_prob = exp_prob * frequency[letter]
            odds = multi_freq[combo] - exp_prob
            Dictionary().add_word(nt_bias, combo, odds)
        return nt_bias

    def poisson_odds(self, sequence, multi_freq, frequency):
        #Calculates nucleotide bias vector using Poisson distribution
        nt_bias = {}
        n = len(sequence)
        for combo in multi_freq.iterkeys():
            exp_prob = 1
            for letter in combo:
                exp_prob = exp_prob * frequency[letter]
            if exp_prob != 0:
                odds = multi_freq[combo] - exp_prob/(math.sqrt(n*exp_prob*(1-exp_prob)))
            else:
                odds = 0
            Dictionary().add_word(nt_bias, combo, odds)
        return nt_bias

    def distance(self, r1, r2):
        #Computes distance between two bias vectors
        sum = 0
        intersect = set(r1.keys()) & set(r2.keys())
        complement = set(r1.keys()) | set(r2.keys()) - intersect
        for key in intersect:
            sum += math.pow((r1[key] - r2[key]),2)
        for key in complement:
            if key in r1:
                sum += math.pow((r1[key]),2)
            elif key in r2:
                sum += math.pow((r2[key]),2)
        distance = math.sqrt(sum)
        return distance

    def odds(self, k, seq1, seq2):
        #Computes nucleotide bias as odds ratio, log of odds ratio, odds difference distribution, and Poisson distribution
        #Distance between bias vectors is computed
        base_freq1 = self.frequency_profile(seq1, 1)[1]
        base_freq2 = self.frequency_profile(seq2, 1)[1]
        freq1 = self.frequency_profile(seq1, k)[1]
        freq2 = self.frequency_profile(seq2, k)[1]
        nt_odds_ratio1 = self.odds_ratio(freq1, base_freq1)
        nt_odds_ratio2 = self.odds_ratio(freq2, base_freq2)
        nt_log_odds1 = self.log_odds_ratio(freq1, base_freq1)
        nt_log_odds2 = self.log_odds_ratio(freq2, base_freq2)
        nt_odds_diff1 = self.odds_difference(freq1, base_freq1)
        nt_odds_diff2 = self.odds_difference(freq2, base_freq2)
        nt_poisson1  = self.poisson_odds(seq1, freq1, base_freq1)
        nt_poisson2 = self.poisson_odds(seq2, freq2, base_freq2)

        nt_or_dist = self.distance(nt_odds_ratio1, nt_odds_ratio2)
        nt_lo_dist = self.distance(nt_log_odds1, nt_log_odds2)
        nt_od_dist = self.distance(nt_odds_diff1, nt_odds_diff2)
        nt_p_dist = self.distance(nt_poisson1, nt_poisson2)
        dist_list = []
        dist_list.append(nt_or_dist)
        dist_list.append(nt_lo_dist)
        dist_list.append(nt_od_dist)
        dist_list.append(nt_p_dist)
        return dist_list

    def cre(self, k, sequence):
        #Computes cumulative relative entropy
        k = int(k)
        profiles, exp = [], {}
        profiles.append(self.frequency_profile(sequence,k))
        profiles.append(self.frequency_profile(sequence,k+1))
        profiles.append(self.frequency_profile(sequence,k+2))
        for word in profiles[2]:
            exp[word] = (profiles[1][word[1:]]*profiles[1][word[:-1]])/profiles[0][word[1:-1]]
        relative_e = self.kldiv(exp, profiles[2])
        return relative_e

class MainFrame(wx.Frame):
    def __init__(self, parent, id, title):
        #Layout of the GUI
        wx.Frame.__init__(self, parent, id, title, size = (845,640))
        panel = wx.Panel(self)
        panel.SetBackgroundColour("#b0c4de")
        sizer = wx.BoxSizer(wx.VERTICAL)
        radio_sizer1 = wx.BoxSizer(wx.VERTICAL)
        spin_sizer = wx.BoxSizer(wx.HORIZONTAL)
        horiz_sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        self.radio_sizer2 = wx.GridBagSizer(3,2)
        
        self.seq_label = wx.StaticText(self, wx.ID_ANY, 'Sequence Type:')
        self.nt_button = wx.RadioButton(self, wx.ID_ANY, 'Nucleotide', style = wx.RB_GROUP)
        self.nt_button.SetValue(True)
        self.aa_button = wx.RadioButton(self, wx.ID_ANY, 'Amino Acid')
        self.gen_button = wx.RadioButton(self, wx.ID_ANY, 'Purine/Pyrimidine')
        radio_sizer1.AddMany([self.seq_label, (self.nt_button, 0, wx.LEFT, 15), (self.gen_button, 0, wx.LEFT, 15), (self.aa_button, 0, wx.LEFT|wx.BOTTOM, 15)])
        
        self.add_button = wx.Button(self, wx.ID_ANY, 'Add Sequences')
        self.add_button.Bind(wx.EVT_BUTTON, self.call_add_sequences)
        self.quick_button = wx.Button(self, wx.ID_ANY, 'Quick Calculation')
        self.quick_button.Bind(wx.EVT_BUTTON, self.call_quick_calc)
        horiz_sizer1.AddMany([(self.add_button, 0, wx.TOP|wx.LEFT|wx.BOTTOM, 5), (self.quick_button, 0, wx.TOP|wx.LEFT, 5), (radio_sizer1, 0, wx.LEFT, 200)])

        self.organismlc = wx.ListCtrl(self, -1, size = (455,564), style = wx.LC_REPORT|wx.LC_HRULES)
        self.organismlc.InsertColumn(0, 'Common Name')
        self.organismlc.InsertColumn(1, 'Accession Number')
        self.organismlc.SetColumnWidth(0,320)
        self.organismlc.SetColumnWidth(1,120)
        
        self.outgroup_label = wx.StaticText(self, wx.ID_ANY, 'Outgroup:')
        self.outgroup = wx.TextCtrl(self, -1, size = (200,20))
        self.outgroup_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.outgroup_sizer.AddMany([(self.outgroup_label), (self.outgroup, 0, wx.LEFT, 5)])
        
        self.tree_label = wx.StaticText(self, wx.ID_ANY, 'Tree Type:')
        self.neighbor_button = wx.RadioButton(self, wx.ID_ANY, 'Neighbor Joining', style = wx.RB_GROUP)
        self.neighbor_button.SetValue(True)
        self.upgma_button = wx.RadioButton(self, wx.ID_ANY, 'UPGMA')
        self.fitch_button = wx.RadioButton(self, wx.ID_ANY, 'Fitch (Additive)')
        self.kitsch_button = wx.RadioButton(self, wx.ID_ANY, 'Kitsch (Ultrametric)')
        self.fitch_button.Bind(wx.EVT_RADIOBUTTON, self.enable_options)
        self.kitsch_button.Bind(wx.EVT_RADIOBUTTON, self.enable_options)
        self.neighbor_button.Bind(wx.EVT_RADIOBUTTON, self.disable_options)
        self.upgma_button.Bind(wx.EVT_RADIOBUTTON, self.disable_options)
        self.fm_button = wx.RadioButton(self, wx.ID_ANY, 'Fitch-Margoliash', style = wx.RB_GROUP)
        self.min_e = wx.RadioButton(self, wx.ID_ANY, 'Minimum Evolution')
        self.fm_button.SetValue(True)
        self.fm_button.Disable()
        self.min_e.Disable()
        
        self.jumble = wx.CheckBox(self, wx.ID_ANY, 'Jumble Sequence Order')
        self.jumble_spin = wx.SpinCtrl(self, style = wx.SP_VERTICAL)
        self.jumble_label = wx.StaticText(self, wx.ID_ANY, 'Number of Jumbles:')
        self.jumble_spin.SetRange(1, 1)
        self.jumble_spin.SetValue(1)
        self.jumble_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.jumble_sizer.AddMany([(self.jumble_label), (self.jumble_spin, 0, wx.TOP, -5)])
        self.power_label = wx.StaticText(self, wx.ID_ANY, 'Power:')
        self.power_spin = wx.SpinCtrl(self, style = wx.SP_VERTICAL)
        self.power_spin.SetRange(0,3)
        self.power_spin.SetValue(2)
        self.power_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.power_sizer.AddMany([self.power_label, (self.power_spin, 0, wx.LEFT, 5)])
        self.g_rearrange = wx.CheckBox(self, wx.ID_ANY, 'Global Rearrangements')
        self.g_rearrange.Disable()
        self.tree_sizer = wx.BoxSizer(wx.VERTICAL)
        self.neighbor_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.fitch_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.neighbor_sizer.AddMany([(self.neighbor_button, 0, wx.LEFT, 10), (self.upgma_button, 0, wx.LEFT, 10)])
        self.fitch_sizer.AddMany([(self.fitch_button, 0, wx.LEFT, 10), (self.kitsch_button, 0, wx.LEFT|wx.RIGHT, 22)])
        self.tree_sizer.AddMany([(self.tree_label, 0, wx.LEFT, 5), self.neighbor_sizer, self.fitch_sizer,
                                 (self.fm_button, 0, wx.LEFT, 35), (self.min_e, 0, wx.LEFT, 35), (self.g_rearrange, 0, wx.LEFT|wx.RIGHT, 35),
                                 (self.power_sizer, 0, wx.LEFT, 35),
                                 (self.jumble, 0, wx.TOP|wx.LEFT, 10), (self.jumble_sizer, 0, wx.TOP|wx.LEFT|wx.BOTTOM, 10)])
        try:
            if os.path.exists('distance.db'):
                connection = sqlite3.connect('distance.db')
                cursor = connection.cursor()
                cursor.execute('''SELECT sequence_name, accession FROM genomes''')
                num_items_org = self.organismlc.GetItemCount()
                for row in sorted(cursor, reverse=True):
                    self.organismlc.InsertStringItem(num_items_org, row[0])
                    self.organismlc.SetStringItem(num_items_org, 1, row[1])
        except sqlite3.OperationalError:
            pass

        self.or_button = wx.RadioButton(self, wx.ID_ANY, 'Odds Ratio', style = wx.RB_GROUP)
        self.lo_button = wx.RadioButton(self, wx.ID_ANY, 'Log Odds')
        self.od_button = wx.RadioButton(self, wx.ID_ANY, 'Odds Difference')
        self.po_button = wx.RadioButton(self, wx.ID_ANY, 'Poisson Odds')
        self.js_button = wx.RadioButton(self, wx.ID_ANY, 'Jensen-Shannon')
        self.mv_button = wx.RadioButton(self, wx.ID_ANY, 'Moment Vector')
        self.nv_button = wx.RadioButton(self, wx.ID_ANY, 'Natural Vector')
        self.or_button.SetValue(True)
        self.or_button.Bind(wx.EVT_RADIOBUTTON, self.select_other)
        self.lo_button.Bind(wx.EVT_RADIOBUTTON, self.select_other)
        self.od_button.Bind(wx.EVT_RADIOBUTTON, self.select_other)
        self.po_button.Bind(wx.EVT_RADIOBUTTON, self.select_other)
        self.js_button.Bind(wx.EVT_RADIOBUTTON, self.select_js)
        self.mv_button.Bind(wx.EVT_RADIOBUTTON, self.select_mv)
        self.nv_button.Bind(wx.EVT_RADIOBUTTON, self.select_nv)
        self.graph_button = wx.Button(self, wx.ID_ANY, 'Graph Sequences')
        self.graph_button.Bind(wx.EVT_BUTTON, self.call_graph)
        self.retrieve_button = wx.Button(self, wx.ID_ANY, 'Retrieve Distances')
        self.retrieve_button.Bind(wx.EVT_BUTTON, self.call_retrieve_distances)
        
        self.graph_sizer = wx.BoxSizer(wx.VERTICAL)
        self.graph_sizer.AddMany([self.graph_button])
        self.label_sizer = wx.BoxSizer(wx.VERTICAL)
        self.calc_label = wx.StaticText(self, wx.ID_ANY, 'Calculation Type:')
        self.freq_label = wx.StaticText(self, wx.ID_ANY, 'Frequency Based:')
        self.vector_label = wx.StaticText(self, wx.ID_ANY, 'Vector Based:')
        self.dimension_label = wx.StaticText(self, wx.ID_ANY, 'Dimension:')
        self.dimension_spin = wx.SpinCtrl(self, style=wx.SP_VERTICAL)
        self.dimension_spin.SetRange(2,25)
        self.dimension_spin.SetValue(2)
        self.vector_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.vector_button_sizer = wx.BoxSizer(wx.VERTICAL)
        self.dimension_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.dimension_sizer.AddMany([self.dimension_label, self.dimension_spin])
        self.vector_button_sizer.AddMany([(self.mv_button, 1, wx.LEFT, 40), (self.nv_button, 1, wx.LEFT, 40), (self.dimension_sizer, 0, wx.LEFT, 70)])
        self.vector_sizer.AddMany([self.vector_button_sizer, (self.graph_sizer,0, wx.LEFT, 5)])
        self.radio_sizer2.AddMany([(self.or_button, (0,0)), (self.lo_button, (0,1)), (self.od_button, (1,0)), (self.po_button, (1,1)), (self.js_button, (2,0))])
        self.label_sizer.AddMany([(self.calc_label, 1, wx.LEFT, 5), (self.freq_label, 1, wx.LEFT, 20),(self.radio_sizer2, 0, wx.LEFT, 40), (self.vector_label, 1, wx.LEFT, 20),
                                  (self.vector_sizer, 0, wx.BOTTOM, 10)])
        
        self.spin = wx.SpinCtrl(self, style=wx.SP_VERTICAL)
        self.spin_label = wx.StaticText(self, wx.ID_ANY, 'Size of k-mers:')
        self.spin.SetRange(2,5)
        self.spin.SetValue(3)
        spin_sizer.AddMany([(self.spin_label, 0, wx.TOP|wx.LEFT, 0), (self.spin, 0, wx.BOTTOM, 20), ])
        
        self.add_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_addseqs)
        self.quick_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_quickcalc)
        self.outgroup.Bind(wx.EVT_ENTER_WINDOW, self.enter_outgroup)
        self.outgroup_label.Bind(wx.EVT_ENTER_WINDOW, self.enter_outgroup)
        self.spin_label.Bind(wx.EVT_ENTER_WINDOW, self.enter_kmer)
        self.spin.Bind(wx.EVT_ENTER_WINDOW, self.enter_kmer)
        self.retrieve_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_retrieve)
        self.nt_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_nucleotide)
        self.gen_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_pp)
        self.aa_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_aminoacid)
        self.jumble.Bind(wx.EVT_ENTER_WINDOW, self.enter_jumble)
        self.jumble_label.Bind(wx.EVT_ENTER_WINDOW, self.enter_jumble_spin)
        self.jumble_spin.Bind(wx.EVT_ENTER_WINDOW, self.enter_jumble_spin)
        self.neighbor_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_neighbor)
        self.upgma_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_upgma)
        self.fitch_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_fitch)
        self.kitsch_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_kitsch)
        self.fm_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_fm)
        self.min_e.Bind(wx.EVT_ENTER_WINDOW, self.enter_mine)
        self.g_rearrange.Bind(wx.EVT_ENTER_WINDOW, self.enter_global)
        self.or_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_oddsratio)
        self.od_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_oddsdiff)
        self.lo_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_logodds)
        self.po_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_poisson)
        self.js_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_js)
        self.mv_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_moment)
        self.nv_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_natural)
        self.dimension_label.Bind(wx.EVT_ENTER_WINDOW, self.enter_dimension)
        self.dimension_spin.Bind(wx.EVT_ENTER_WINDOW, self.enter_dimension)
        self.graph_button.Bind(wx.EVT_ENTER_WINDOW, self.enter_graph)
        self.power_label.Bind(wx.EVT_ENTER_WINDOW, self.enter_power)
        self.power_spin.Bind(wx.EVT_ENTER_WINDOW, self.enter_power)
        
        sizer.AddMany([horiz_sizer1,(self.label_sizer, 0, wx.LEFT, 470), (spin_sizer, 1, wx.LEFT|wx.RIGHT, 475),
                       (self.tree_sizer, 0, wx.LEFT, 470),(self.outgroup_sizer, 1, wx.LEFT, 480), 
                       (self.retrieve_button, 0, wx.LEFT, 475), (self.organismlc, 0, wx.TOP, -560)])
        
        self.sb = self.CreateStatusBar()
        self.SetMinSize(self.GetSize())
        self.SetMaxSize(self.GetSize())
        
        panel.SetSizerAndFit(sizer)
        self.Centre()
        self.Show(True)
    
    #Below methods show information about an item in the status bar when the user hovers over that item 
    def enter_addseqs(self, event):
        self.sb.SetStatusText("Add sequences to database from FASTA file")
        event.Skip()
    def enter_quickcalc(self, event):
        self.sb.SetStatusText("Calculate distance for selected calculation only for sequences in FASTA file.  Outgroup should be first sequence in file.")
        event.Skip()
    def enter_outgroup(self, event):
        self.sb.SetStatusText("Enter the accession number of the virus serving as the outgroup")
        event.Skip()
    def enter_kmer(self, event):
        self.sb.SetStatusText("Select the size of the k-mer used in the calculation")
        event.Skip()
    def enter_retrieve(self, event):
        self.sb.SetStatusText("Click to create distance matrix and tree for all selected viruses")
        event.Skip() 
    def enter_nucleotide(self, event):
        self.sb.SetStatusText("Select kmers of nucleotides")
        event.Skip()
    def enter_aminoacid(self, event):
        self.sb.SetStatusText("Select kmers of amino acids")
        event.Skip()
    def enter_pp(self, event):
        self.sb.SetStatusText("Generalize nucleotide sequence to purines and pyrimidines: A & G to 'R' and C & T to 'Y' to potentially reduce strain specific differences")
        event.Skip()
    def enter_jumble(self, event):
        self.sb.SetStatusText("Randomize order of sequences input to tree drawing program")
        event.Skip()
    def enter_jumble_spin(self, event):
        self.sb.SetStatusText("Select number of times to randomize.  Best tree will be selected")
        event.Skip()
    def enter_neighbor(self, event):
        self.sb.SetStatusText("NJ: generates unrooted tree by successive clustering of lineages, setting branch lengths as lineages join. No rearrangement, no evolutionary clock")
        event.Skip()
    def enter_upgma(self, event):
        self.sb.SetStatusText("UPGMA: constructs rooted tree by successive (agglomerative) clustering using average-linkage method of clustering")
        event.Skip()
    def enter_fitch(self, event):
        self.sb.SetStatusText("Fitch: Contructs phylogeny under 'additive tree model' where distances are expected to equal sums of branch lengths between the species")
        event.Skip()  
    def enter_kitsch(self, event):
        self.sb.SetStatusText("Kitsch: Uses the 'ultrametric model' which is the same as the additive model of Fitch but an evolutionary clock is assumed")
        event.Skip()  
    def enter_fm(self, event):
        self.sb.SetStatusText("Uses Fitch-Margolish criterion to fit branch lengths to each topology based on goodness of fit sum of squares")
        event.Skip()
    def enter_mine(self,event):
        self.sb.SetStatusText("Minimum Evolution: Uses Fitch-Margoliash criterion to fit branch lengths to each topology, but chooses topologies based on total branch length")
        event.Skip()
    def enter_global(self,event):
        self.sb.SetStatusText("Global rearrangements: after last species added, each possible group is removed and re-added to improve accuracy.  ~3X Runtime")
        event.Skip()
    def enter_oddsratio(self, event):
        self.sb.SetStatusText("Calculate k-mer bias of an organism as an odds ratio.  The distance is Euclidean distance of the bias of two organisms.")
        event.Skip()
    def enter_oddsdiff(self, event):
        self.sb.SetStatusText("Calculate k-mer bias of an organism as an odds difference.  Similar to Poisson without the normalization.  The Euclidean distance is calculated.")
        event.Skip()
    def enter_logodds(self, event):
        self.sb.SetStatusText("Calculate k-mer bias of an organism as the logarithm of the odds ratio.  The distance is Euclidean distance of the bias of two organisms.")
        event.Skip()
    def enter_poisson(self, event):
        self.sb.SetStatusText("Calculate k-mer bias as the normalized deviation from a Poisson distribution.  The distance is Euclidean distance of the bias of two organisms.")
        event.Skip()
    def enter_js(self, event):
        self.sb.SetStatusText("Calculates Jensen-Shannon distance using Kullback-Leibler divergence from k-mer frequency distribution")
        event.Skip()
    def enter_moment(self, event):
        self.sb.SetStatusText("Construct graphical representation of sequence and characterizes curve using moment vectors.  The Euclidean distance is calculated.")
        event.Skip()
    def enter_natural(self, event):
        self.sb.SetStatusText("Similar to moment vector analysis, but uses natural vectors, which have associated parameters describing the numbers and distribution of k-mers.")
        event.Skip()
    def enter_dimension(self, event):
        self.sb.SetStatusText("Choose the dimension of the vector to use for calculations")
        event.Skip()
    def enter_graph(self, event):
        self.sb.SetStatusText("Graphical representation of DNA/amino acid sequences")
        event.Skip()
    def enter_power(self, event):
        self.sb.SetStatusText("Select power for Fitch and Kitsch. 2 is the default for Fitch-Margoliash.  The power is 0 for UPGMA and Neighbor (implicitly).")
        event.Skip()
        
    def seed_number(self):
        #Generates random odd number of the type 4n+1 to use as a seed for jumbling the sequence order
        num = 0
        while num%4 != 1:
            num = random.randrange(1,2147483647,2)
        return num
              
    def enable_options(self, event):
        #Enables some tree type options
        self.fm_button.Enable()
        self.min_e.Enable()
        if self.fitch_button.GetValue():
            self.g_rearrange.Enable()
        else:
            self.g_rearrange.Disable()
        self.jumble_spin.SetValue(10)
        self.jumble_spin.SetRange(1,50)  
    def disable_options(self, event):
        #Disables some tree type options
        self.fm_button.Disable()
        self.min_e.Disable()
        self.g_rearrange.Disable()
        if self.neighbor_button.GetValue():
            self.jumble_spin.SetValue(1)
            self.jumble_spin.SetRange(1,1)
        else:
            self.jumble_spin.SetValue(10)
            self.jumble_spin.SetRange(1,50)

    def select_js(self, event):
        #Set appropriate ranges when Jensen-Shannon radio button selected
        if not self.gen_button.Enabled:
            self.gen_button.Enable()
        self.spin.SetRange(3,8)
        self.spin.SetValue(3)
        self.Update()  
    def select_mv(self, event):
        #Disable certain options and set appropriate ranges when moment vector radio button selected
        self.nt_button.SetValue(True)
        self.dimension_spin.SetRange(2,50)
        if self.gen_button.Enabled:
            self.gen_button.Disable()
        self.Update()  
    def select_nv(self, event):
        #Disable certain options and set appropriate ranges when natural vector radio button selected
        self.nt_button.SetValue(True)
        self.dimension_spin.SetRange(2,25)
        if self.gen_button.Enabled:
            self.gen_button.Disable()
        self.spin.SetRange(1,5)
        self.spin.SetValue(1)
        self.Update()
    def select_other(self, event):
        #Reenable the generalized button if vector buttons are not selected
        if not self.gen_button.Enabled:
            self.gen_button.Enable()
        self.spin.SetRange(2,5)
        self.spin.SetValue(3)
        self.Update()
        
    def get_selected_items(self,list_control):
        #Get selected items in list control
        selection = []
        current = -1
        while True:
            next = self.getNextSelected(list_control, current)
            if next == -1:
                return selection
            selection.append(next)
            current = next
    def getNextSelected(self,list_control, current):
        #Get next selected item in list control
        return list_control.GetNextItem(current,wx.LIST_NEXT_ALL,wx.LIST_STATE_SELECTED)
    
    def call_add_sequences(self, event):
        #Call add sequences method from input class input sequences to database
        I = Input()
        file = I.add_sequences()
        if file != "":
            reads = I.read_seqs_from_file(file)
            for seq_name, accession in itertools.izip(reads[0],reads[3]):
                num = 0
                org_name = ""
                for org in range(self.organismlc.GetItemCount()): #@UnusedVariable
                    org_name = self.organismlc.GetItemText(num)
                    if seq_name > org_name:
                        num += 1
                if seq_name != org_name:
                    self.organismlc.InsertStringItem(num, seq_name)
                    self.organismlc.SetStringItem(num, 1, accession)
        self.Update()
    def call_quick_calc(self, event):
        #Calls the Quick Calc method
        I = Input()
        file = I.input_file()
        if file != "":
            reads = I.read_seqs_from_file(file)
        accession_dict = {}
        try:
            for x in xrange(len(reads)):
                accession_dict[reads[3][x]] = reads[0][x]
            if self.nt_button.GetValue():
                type = 'nucleotide'
            elif self.aa_button.GetValue():
                type = 'aminoacid'
            elif self.gen_button.GetValue():
                type = 'generalized'  
            if self.or_button.GetValue():
                calculation = 'oddsratio'
            elif self.lo_button.GetValue():
                calculation = 'logodds'
            elif self.od_button.GetValue():
                calculation = 'oddsdiff'
            elif self.po_button.GetValue():
                calculation = 'poisson'
            elif self.js_button.GetValue():
                calculation = 'jensen_shannon'
            elif self.mv_button.GetValue():
                calculation = 'moment'
            elif self.nv_button.GetValue():
                calculation = 'natural'
            k = self.spin.GetValue()
            Database().quick_calc(int(k), reads, type, calculation)
            self.draw_tree(accession_dict)
        except UnboundLocalError:
            pass
    def call_graph(self, event):
        #Calls graph method from database class for frequency calculations and draw method from vector class in genome space module for vector calculations
        tmp_list, org_list, org_dict = [], [], {}
        num = 0
        list = self.get_selected_items(self.organismlc)
        connection = sqlite3.connect('distance.db')
        cursor = connection.cursor()
        cursor.execute('''SELECT sequence_name, accession, nt_sequence FROM genomes''')
        for row in cursor:
            tmp_list.append(str(row[0] + '\t' + row[1] + '\t' + row[2]))
        tmp_list = sorted(tmp_list)
        for row in tmp_list:
            if num in list:
                org_list.append(row.split('\t')[1])
                org_dict[row.split('\t')[0]] = row.split('\t')[2]
            num += 1
        V = gs.Vectors()
        if not self.mv_button.GetValue() and not self.nv_button.GetValue():
            if self.nt_button.GetValue():
                type = 'nucleotide'
            elif self.aa_button.GetValue():
                type = 'aminoacid'
            elif self.gen_button.GetValue():
                type = 'generalized'
                
            if self.or_button.GetValue():
                calculation = 'oddsratio'
            elif self.lo_button.GetValue():
                calculation = 'logodds'
            elif self.od_button.GetValue():
                calculation = 'oddsdiff'
            elif self.po_button.GetValue():
                calculation = 'poisson'
            elif self.js_button.GetValue():
                calculation = 'jensen_shannon'
            DB = Database()
            dist = DB.retrieve_distances(self.spin.GetValue(), org_list, type, calculation, False)
            DB.graph(dist[0], dist[1])
        else:
            if self.nt_button.GetValue():
                type = 'nt'
            elif self.gen_button.GetValue():
                type = 'nt'
            elif self.aa_button.GetValue():
                type = 'aa'
            if self.mv_button.GetValue():
                calculation = 'moment'
            elif self.nv_button.GetValue():
                calculation = 'natural'
            dimension = int(self.dimension_spin.GetValue())
            k = int(self.spin.GetValue())
            V.draw(org_dict, type, calculation, dimension, k)
    def call_retrieve_distances(self, event):
        #Calls retrieve distances method from Database class
        try:
            org_list, tmp_list, alist = [], [], []
            num = 0
            accession_dict = {}
            list = self.get_selected_items(self.organismlc)
            connection = sqlite3.connect('distance.db')
            cursor = connection.cursor()
            cursor.execute('''SELECT sequence_name, accession FROM genomes''')
            for row in cursor:
                tmp_list.append(str(row[0] + '\t' + row[1]))
            for x in sorted(tmp_list):
                alist.append(x)
            for row in alist:
                if num in list:
                    org_list.append(row.split('\t')[1])
                    accession_dict[row.split('\t')[1]] = row.split('\t')[0]
                num += 1
            if self.nt_button.GetValue():
                type = 'nucleotide'
            elif self.aa_button.GetValue():
                type = 'aminoacid'
            elif self.gen_button.GetValue():
                type = 'generalized'
                
            if self.or_button.GetValue():
                calculation = 'oddsratio'
            elif self.lo_button.GetValue():
                calculation = 'logodds'
            elif self.od_button.GetValue():
                calculation = 'oddsdiff'
            elif self.po_button.GetValue():
                calculation = 'poisson'
            elif self.js_button.GetValue():
                calculation = 'jensen_shannon'
            elif self.mv_button.GetValue():
                calculation = 'moment'
            elif self.nv_button.GetValue():
                calculation = 'natural'
            Database().retrieve_distances(self.spin.GetValue(), org_list, type, calculation, True)
            self.draw_tree(accession_dict)
        except sqlite3.OperationalError:
            pass
        
    def draw_tree(self, accession_dict):
        #Use Phylip distance matrix programs: Neighbor, Fitch, or Kitsch to draw various types of phylogenetic trees
        out = 1
        string = ""
        #Identify number of outgroup
        if self.outgroup.GetValue() != "":
            og = str(self.outgroup.GetValue())
            num = 0
            with open ('infile') as file:
                for line in file:
                    if line[:9] == og[:9]:
                        out = num
                    num += 1
        #Provide options to Phylip distance matrix programs
        if self.neighbor_button.GetValue():
            if self.outgroup.GetValue():
                string += 'O\n' + str(out) + '\n'
            if self.jumble.GetValue():
                string += 'J\n' + str(self.seed_number()) + '\n' + str(self.jumble_spin.GetValue()) + '\n'
            string += 'Y\nR\n'
            f = open('options.txt', 'w')
            f.write(string)
            f.close()
            os.system("./neighbor < options.txt")
        elif self.upgma_button.GetValue():
            string += "N\n"
            if self.outgroup.GetValue():
                string += 'O\n' + str(out) + '\n'
            if self.jumble.GetValue():
                string += 'J\n' + str(self.seed_number()) + '\n' + str(self.jumble_spin.GetValue()) + '\n'
            string += 'Y\nR\n'
            f = open('options.txt', 'w')
            f.write(string)
            f.close()
            os.system("./neighbor < options.txt")
        elif self.fitch_button.GetValue():
            if self.min_e.GetValue():
                string += 'D\n'
            if self.outgroup.GetValue():
                string += 'O\n' + str(out) + '\n'
            if self.g_rearrange.GetValue():
                string += 'G\n'
            string += 'P\n' + str(self.power_spin.GetValue()) + '\n'
            if self.jumble.GetValue():
                string += 'J\n' + str(self.seed_number()) + '\n' + str(self.jumble_spin.GetValue()) + '\n'
            string += 'Y\nR\n'
            f = open('options.txt', 'w')
            f.write(string)
            f.close()
            os.system("./fitch < options.txt")
        elif self.kitsch_button.GetValue():
            if self.min_e.GetValue():
                string += 'D\n'
            string += 'P\n' + str(self.power_spin.GetValue()) + '\n'
            if self.jumble.GetValue():
                string += 'J\n' + str(self.seed_number()) + '\n' + str(self.jumble_spin.GetValue()) + '\n'
            string += 'Y\nR\n'
            f = open('options.txt', 'w')
            f.write(string)
            f.close()
            os.system("./kitsch < options.txt")
        os.remove('infile')
        os.remove('outfile')
        os.remove('options.txt')
        lines = []
        with open('outtree') as f:
            for line in f:
                string = line
                for accession in accession_dict.iterkeys():
                    for m in re.finditer(accession.split('.')[0], string):
                        tmp1 = string[:m.start()]
                        tmp2 = string[m.end():]
                        insert = accession_dict[accession].strip(' ')
                        string = tmp1 + insert + tmp2
                lines.append(string)
        os.remove('outtree')
        file = open('outtree.newick', 'w')
        for line in lines:
            file.write(line)
        file.close()
        t = ThreadClass()
        t.start()
            
class ThreadClass(threading.Thread):
    #Create thread to launch Archaeopteryx for drawing the tree
    def run(self):
        os.system("java -cp forester.jar org.forester.archaeopteryx.Archaeopteryx -c _aptx_configuration_file outtree.newick")
        try:
            os.remove('outtree.newick')
        except OSError:
            pass
        
class FastMap: 
    def __init__(self, dist,verbose=False): 
        if dist.max()>1:
            dist/=dist.max()
        self.dist=dist
        self.verbose=verbose

    def _furthest(self, o): 
        mx=-1000000
        idx=-1
        for i in range(len(self.dist)): 
            d=self._dist(i,o, self.col)
            if d>mx: 
                mx=d
                idx=i
        return idx

    def _pickPivot(self):
        #Find the two most distant points
        o1=random.randint(0, len(self.dist)-1)
        o2=-1
        i=DISTANCE_ITERATIONS
        while i>0: 
            o=self._furthest(o1)
            if o==o2: break
            o2=o
            o=self._furthest(o2)
            if o==o1: break
            o1=o
            i-=1
        self.pivots[self.col]=(o1,o2)
        return (o1,o2)

    def _map(self, K): 
        if K==0: return 
        px,py=self._pickPivot()
        if self.verbose: pass
        if self._dist(px,py,self.col)==0: 
            return 
        for i in range(len(self.dist)):
            self.res[i][self.col]=self._x(i, px,py)
        self.col+=1
        self._map(K-1)

    def _x(self,i,x,y):
        #Project the i'th point onto the line defined by x and y
        dix=self._dist(i,x,self.col)
        diy=self._dist(i,y,self.col)
        dxy=self._dist(x,y,self.col)
        return (dix + dxy - diy) / 2*math.sqrt(dxy)

    def _dist(self, x,y,k): 
        #Recursively compute the distance based on previous projections
        if k==0: return self.dist[x,y]**2
        rec=self._dist(x,y, k-1)
        resd=(self.res[x][k] - self.res[y][k])**2
        return rec-resd

    def map(self, K): 
        self.col=0
        self.res=scipy.zeros((len(self.dist),K))
        self.pivots=scipy.zeros((K,2),"i")
        self._map(K)
        return self.res

def main():
    app = wx.App(redirect=False)
    MainFrame(None, -1, 'Virus Classifier')
    app.MainLoop()   
        
if __name__ == "__main__":
    main(*sys.argv[1:])