# HW3-python

## Kmer spectrum creater and visualizer  

class Kmer_spectrum  
_Input parameters:_ kmer size (k, int), quality threshold (q, int)  
_Methods:_   
* create_spectrum (parameter: filename.fastq)
* calc_frequency: calculate the abundance of frequencies across dictionary, created by _create_spectrum_
* visualise_spectrum (parameter: list of four integers [xmin,xmax,ymin,ymax] - scales of axes in spectrum visualisation): plot frequncies. helps to decide the sequncing errors / valuable data border
* calculate_genome_size (parameter: noize threshold): calculate approximate genome size based on kmers frequencies.

## Results for __test_kmer.fastq.__
```
test = Kmer_spectrum(15)  
test.create_spectrum('test_kmer.fastq')  
test.calc_frequency()  
test.visualize_spectrum([0,250, 1000, 1000000])  
test.calculate_genome_size(50)  ```
