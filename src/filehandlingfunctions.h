int nColInFile(const char *fname);
size_t readFrequencyFile(const char *fname,std::vector<double> &ret);
bgl findOverlap(const char* glfileA, const char* glfileB, const char* freqfile, std::vector<double> &freq);
bgl readBeagle(const char* fname);
gzFile openFileGz(const char* a,const char* b,const char *mode,std::vector<char *>& dumpedFiles);
FILE *openFile(const char* a,const char* b,std::vector<char *>& dumpedFiles);
