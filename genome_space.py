'''
Created on Jun 22, 2011
@author: Douglas Crandell
'''
from __future__ import division
from Bio.Seq import Seq
import matplotlib #@UnresolvedImport
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt #@UnresolvedImport
import matplotlib.font_manager #@UnresolvedImport
import math
import itertools
        
class Vectors():
    def nt_graph(self, sequence):
        #Construct a DNA sequence graph in two quadrants of Cartesian coordinate system with pyrimidines (C and T)
        #in the first quadrant and purines (A and G) in thr fourth quadrant
        y = 0
        yvalues = []
        yvalues.append(y)
        for letter in sequence:
            if letter == 'C':
                y += 2/3
            elif letter == 'T':
                y += 1/3
            elif letter == 'A':
                y -= 1/3
            elif letter == 'G':
                y -=2/3
            else:
                pass
            yvalues.append(y)
        return yvalues
    
    def aa_graph(self, sequence):
        #Construct an amino acid sequence graph on two quadrants of the Cartesian coordinate system.  Y-coordinate values
        #of 20 amino acids is based on amino acid hydrophobicity scale values.
        y = 0
        yvalues = []
        yvalues.append(y)
        sequence = str(Seq(sequence).translate())
        for letter in sequence:
            if letter == 'W':
                y += 12/13
            elif letter == 'I':
                y += 11/13
            elif letter == 'F':
                y += 10/13
            elif letter == 'L':
                y += 9/13
            elif letter == 'C':
                y += 8/13
            elif letter == 'M':
                y += 7/13
            elif letter == 'V':
                y += 6/13
            elif letter == 'Y':
                y += 5/13
            elif letter == 'P':
                y += 4/13
            elif letter == 'A':
                y += 3/13
            elif letter == 'T':
                y += 2/13
            elif letter == 'H':
                y += 1/13
            elif letter == 'G':
                y += 0
            elif letter == 'S':
                y -= 1/8
            elif letter == 'Q':
                y -= 2/8
            elif letter == 'N':
                y -= 3/8
            elif letter == 'E':
                y -= 4/8
            elif letter == 'D':
                y -= 5/8
            elif letter == 'K':
                y -= 6/8
            elif letter == 'R':
                y -= 7/8
            else:
                pass
            yvalues.append(y)
        return yvalues
                
    def moment_vector(self, values):
        #Construct moment vector from sequence graph vector.
        n = len(values)
        mv = []
        for j in range(1, 51):
            moment = 0
            for x in xrange(1, n+1):
                moment += math.pow((x- values[x-1]),j)/(math.pow(n,j))
            mv.append(moment)
        return mv

    def frequency_profile(self, sequence, k):
        #Get frequency distribution for nucleotide/amino acid kmers
        k = int(k)
        count = 0
        string = ""
        profile = {}
        for letter in sequence:
            count += 1
            if count < k:
                string = string + letter
            elif count == k:
                string = string + letter
                profile[string] = profile.get(string, []) + [count-k]
            else:
                string = string[1:] + letter
                profile[string] = profile.get(string, []) + [count-k]
        return profile
    
    def distance_mean(self, profile):
        #Compute average distance from origin for kmers in frequency distribution
        means = {}
        for pattern in profile.iterkeys():
            means[pattern] = sum(profile[pattern])/len(profile[pattern])
        return means
    
    def natural_vector(self, sequence, profile):
        #Compute natural vector from normalized central moments to generate numerical characterization of sequence
        nv = []
        means = self.distance_mean(profile)
        for pattern in profile.iterkeys():
            central_moment = 0
            length = 0
            if len(profile[pattern]) > 25:
                length = 26
            else:
                length = len(profile[pattern]) + 1
            for j in xrange(1, length):
                if j == 1:
                    central_moment = 0
                else:
                    central_moment += math.pow(sum(profile[pattern])-means[pattern], j)/(math.pow(len(profile[pattern]),j-1)*math.pow(len(sequence),j-1))
                nv.append(central_moment)
        return nv
    
    def vector_distance(self, v1, v2):
        #Compute the distance between two vectors
        distance = 0
        for moment1, moment2 in itertools.izip(v1,v2):
            distance += math.pow((float(moment1) - float(moment2)),2)
        distance = math.sqrt(distance)
        return distance
    
    def draw(self, dict, type, calculation, dimension, k):
        #Plots two graphs.  The DNA or amino acid graphical representation and the moment or natural vector.
        #The first component of the moment or natural vector will be plotted on the x-axis and the selected
        #dimension will be plotted on the y-axis.
        prop = matplotlib.font_manager.FontProperties(size=8)
        font = {'fontsize':12}
        labelfont = {'fontsize':9}
        colors = 'bgrcmyk'
        ax = plt.subplot(211)
        num = 0
        x, y, annotes = [], [], []
        if type == 'nt':
            base = 'nucleotide'
            for name, seq in dict.items():
                yvalues = self.nt_graph(seq)
                if num == len(colors):
                    num = 0
                plt.plot(yvalues, c=colors[num], label=name)
                plt.title('Nucleotide Composition', font)
                plt.xlabel('Number of Nucleotides', labelfont)
                plt.ylabel('Cumulative y-value of composition vector',labelfont)
                num += 1
                if calculation == 'moment':
                    mv = self.moment_vector(yvalues)
                    x.append(mv[0])
                    y.append(mv[dimension-1])
                    annotes.append(name)
                else:
                    fp = self.frequency_profile(seq, k)
                    nv = self.natural_vector(seq, fp)
                    x.append(nv[1])
                    y.append(nv[dimension-1])
                    annotes.append(name)
        else:
            base = 'amino acid'
            for name, seq in dict.items():
                yvalues = self.aa_graph(seq)
                if num == len(colors):
                    num = 0
                plt.plot(yvalues, c=colors[num], label=name)
                plt.title('Amino Acid Composition', font)
                plt.xlabel('Number of Amino Acids', labelfont)
                plt.ylabel('Cumulative y-value of composition vector', labelfont)
                num += 1
                if calculation == 'moment':
                    mv = self.moment_vector(yvalues)
                    x.append(mv[0])
                    y.append(mv[dimension-1])
                    annotes.append(name)
                else:
                    fp = self.frequency_profile(seq, k)
                    nv = self.natural_vector(seq, fp)
                    x.append(nv[1])
                    y.append(nv[dimension-1])
                    annotes.append(name)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='best', bbox_to_anchor=(0.5, 1.05),
          ncol=1, prop=prop)
        plt.subplot(212)
        plt.scatter(x,y, c=colors)
        if calculation == 'moment':
            plt.title('Moment Vector: %s, dimension=%s' %(base,dimension), font)
            plt.xlabel('M1', labelfont)
            plt.ylabel('M%s' %str(dimension), labelfont)
        else:
            plt.title('Natural Vector: %s, dimension=%s' %(base,dimension), font)
            plt.xlabel('N1', labelfont)
            plt.ylabel('N%s' %str(dimension), labelfont)
        plt.subplots_adjust(hspace = 0.5)
        af = AnnoteFinder(x,y,annotes)
        plt.connect('button_press_event', af)
        plt.show()
        
class AnnoteFinder():
    def __init__(self, xdata, ydata, annotes, axis=None, xtol=None, ytol=None):
        self.data = zip(xdata, ydata, annotes)
        if xtol is None:
            xtol = ((max(xdata) - min(xdata))/float(len(xdata)))/2
        if ytol is None:
            ytol = ((max(ydata) - min(ydata))/float(len(ydata)))/2
        self.xtol = xtol
        self.ytol = ytol
        if axis is None:
            self.axis = plt.gca()
        else:
            self.axis= axis
        self.drawnAnnotations = {}
        self.links = []

    def distance(self, x1, x2, y1, y2):
        #Return the distance between two points
        return(math.sqrt( (x1 - x2)**2 + (y1 - y2)**2 ))
    
    def __call__(self, event):
        #Method called when button press event is registered on plot
        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata
            if (self.axis is None) or (self.axis==event.inaxes):
                annotes = []
                for x,y,a in self.data:
                    if  (clickX-self.xtol < x < clickX+self.xtol) and  (clickY-self.ytol < y < clickY+self.ytol) :
                        annotes.append((self.distance(x,clickX,y,clickY),x,y, a) )
                if annotes:
                    annotes.sort()
                    distance, x, y, annote = annotes[0] #@UnusedVariable
                    self.drawAnnote(event.inaxes, x, y, annote)
                    for l in self.links:
                        l.drawSpecificAnnote(annote)
        
    def drawAnnote(self, axis, x, y, annote):
        #Draw the annotation on the plot
        if self.drawnAnnotations.has_key((x,y)):
            markers = self.drawnAnnotations[(x,y)]
            for m in markers:
                m.set_visible(not m.get_visible())
            self.axis.figure.canvas.draw()
        else:
            t = axis.text(x,y, " - %s"%(annote))
            m = axis.scatter([x],[y], marker='d', c='r', zorder=100)
            self.drawnAnnotations[(x,y)] =(t,m)
            self.axis.figure.canvas.draw()
    
    def drawSpecificAnnote(self, annote):
        annotesToDraw = [(x,y,a) for x,y,a in self.data if a==annote]
        for x,y,a in annotesToDraw:
            self.drawAnnote(self.axis, x, y, a)