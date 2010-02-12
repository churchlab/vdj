enum nucleotide {A,T,C,G};

unsigned short getNucleotide(char ch){
       switch(ch) {
        case 'a':
        case 'A': return A;
        case 't':
        case 'T': return T;
        case 'c':
        case 'C': return C;
        case 'g': 
        case 'G': return G;
        default:  return 5;
    }
}
