#include "Functions/functions.h"

int main(){
    std::string filename = "Data/test/this/directory/somefile.dat";

    GCT::create_directory_to_file(filename);

    return 0;
}