from abc import ABC, abstractmethod

# 1. Data Abstraction & Polymorphism
class DNASequence(ABC):
    # Abstract class to define basic DNA sequence operations
    def __init__(self, sequence):
        self._sequence = sequence  # Protected attribute for encapsulation

    @abstractmethod
    def get_type(self):
        pass  # Abstract method for identifying sequence type

    def get_sequence(self):
        # Public method to access sequence
        return self._sequence

    @abstractmethod
    def calculate_gc_content(self):
        # Abstract method for polymorphic GC content calculation
        pass

# 2. Inheritance & Polymorphism
class CodingDNA(DNASequence):
    def __init__(self, sequence, gene_name):
        super().__init__(sequence)
        self.gene_name = gene_name

    def get_type(self):
        return "Coding DNA"

    def calculate_gc_content(self):
        # Calculate GC content specific to coding DNA
        gc_content = (self._sequence.count("G") + self._sequence.count("C")) / len(self._sequence) * 100
        return f"GC Content for {self.gene_name}: {gc_content:.2f}%"

class NonCodingDNA(DNASequence):
    def get_type(self):
        return "Non-Coding DNA"

    def calculate_gc_content(self):
        # Calculate GC content for non-coding regions
        gc_content = (self._sequence.count("G") + self._sequence.count("C")) / len(self._sequence) * 100
        return f"GC Content for non-coding region: {gc_content:.2f}%"

# 3. Encapsulation
class DNALibrary:
    def __init__(self, library_name):
        self.__library_name = library_name  # Private attribute for encapsulation
        self.__sequences = []  # Private list to store DNA sequences

    def add_sequence(self, dna_sequence):
        # Public method to add a DNA sequence to the library
        self.__sequences.append(dna_sequence)

    def get_sequences_info(self):
        # Public method to retrieve information about all sequences
        return [(seq.get_type(), seq.get_sequence()) for seq in self.__sequences]

    def get_gc_contents(self):
        # Public method to retrieve GC content for each sequence
        return [seq.calculate_gc_content() for seq in self.__sequences]

    def get_library_name(self):
        # Public method to access library name
        return self.__library_name

# Testing the functionality
def main():
    # Creating instances of DNA sequences with polymorphism
    coding_dna = CodingDNA("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", "GeneA")
    non_coding_dna = NonCodingDNA("GCGTGCATGTCAGGCTCGTAGCT")

    # Creating a DNA library and adding sequences to it
    dna_library = DNALibrary("Genomics Lab Library")
    dna_library.add_sequence(coding_dna)
    dna_library.add_sequence(non_coding_dna)

    # Displaying library name using encapsulation
    print(f"Library Name: {dna_library.get_library_name()}\n")

    # Displaying information about each DNA sequence in the library
    print("DNA Sequences in the Library:")
    for seq_type, seq in dna_library.get_sequences_info():
        print(f"Type: {seq_type}, Sequence: {seq}")

    # Calculating and displaying GC content for each DNA sequence
    print("\nGC Content for DNA Sequences:")
    for gc_content in dna_library.get_gc_contents():
        print(gc_content)

# Run the main function
if __name__ == "__main__":
    main()
