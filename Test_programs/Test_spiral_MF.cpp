#include "Bfield/class_bfield.h"

#include <fstream>
#include <vector>

int main(void){

    Bfield bfield;
    std::array<double, 3> pos, pos2;

    std::ofstream file, file2;

    file.open("Data/spiral_MF.dat");
    file2.open("Data/spiral_MF2.dat");

    for(double x = -500; x <= 500; x+=50){
        for(double y = -500; y <= 500; y+=50){
            
            pos = { x, y, 0 };
            pos2 = { x+25, y+25, 0 };

            bfield.generate_bfield_at_point(0.0, pos);
            file << pos[0] << ' ' << pos[1] << ' ' << bfield.B[0] << ' ' << bfield.B[1] << '\n';
            
            bfield.generate_bfield_at_point(0.0, pos2);
            file2 << pos2[0] << ' ' << pos2[1] << ' ' << bfield.B[0] << ' ' << bfield.B[1] << '\n';
        }
    }

    file.close();
    file.close();

    std::cout << "Wrote spiral field to Data/spiral_MF.dat \n";

}