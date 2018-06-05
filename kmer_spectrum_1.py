
# coding: utf-8

# In[16]:


from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

class Kmer_spectrum:
    def __init__ (self,k,q=1):
        self.kmer_size = k
        self.quality_thres = q
        self.spectrum = {}
        
        
    def create_spectrum(self, filename):
        spectrum = {}
        with open (filename, 'r') as input:
            for record in SeqIO.parse(input, 'fastq'):
                read_length = len(record.seq)
                for i in range(0,read_length-self.kmer_size):
                    current_kmer = str(record.seq[i:i+self.kmer_size])
                    read_quality = record.letter_annotations['phred_quality'][i:i+self.kmer_size]
                    if not len(list(filter(lambda x: x <= self.quality_thres, read_quality))): 
                        if current_kmer not in spectrum:
                            spectrum.update({current_kmer:1})
                        else:
                            spectrum[current_kmer] +=  1
        self.spectrum = spectrum
                            
    def calc_frequency(self):
        frequencies = list(self.spectrum.values())
        final_spectrum = {}
        for i in frequencies:
            if i in final_spectrum:
                final_spectrum[i] += 1
            else:
                final_spectrum.update({i:1})
        self.final_spectrum = final_spectrum
        self.frequencies = frequencies
        #print(final_spectrum)
    
    def visualize_spectrum(self, axis):
        frequencies2 = []
        for key, value in self.final_spectrum.items():
            frequencies2.append((key,value))
        frequencies2.sort(key=lambda tup: tup[0])
        self.sorted_frequencies = frequencies2
        x,y = zip(*frequencies2)
        x = list(x)
        y =list(y)
        y_log = np.log10(y)
        fig, ax = plt.subplots(figsize=(20,10))
        ax.plot(x,y)
        plt.axis([0,250, 1000, 1000000])
        plt.yscale('log')
        plt.title('Kmer spectrum visualisation')
        plt.xlabel('Frequencies')
        plt.ylabel('Number of Kmers')
        return ax
        
    def calculate_genome_size(self,thres):
        not_trash = [x for x in self.sorted_frequencies if x[0] > thres]
        genome_size = sum([x[0]*x[1] for x in not_trash])/self.kmer_size
        print('Genome size is {:.0f} base pairs'.format(genome_size))
        
get_ipython().run_line_magic('matplotlib', 'inline')


# In[18]:


a = Kmer_spectrum(15)
a.create_spectrum('test_kmer.fastq')
a.calc_frequency()
a.visualize_spectrum([0,250, 1000, 1000000])
a.calculate_genome_size(50)


# In[15]:




