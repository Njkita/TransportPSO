#pragma once

#include <iostream>
#include <fstream>

const int quantity = 15;
const int mine = 10;
const int terminal = 5;

class Graph
{
public:
    double roads[quantity][quantity]{},
        wantMine[mine]{},
        wantTerminal[terminal]{};

    Graph(const char* filename) 
    {
        FILE* file;
        errno_t returnValue = fopen_s(&file, filename, "r");

        std::ifstream in;
        in.open(filename);

        if (in.is_open())
        {
            for (int i = 0; i < quantity; i++)
                for (int j = 0; j < quantity; j++)
                    in >> roads[i][j];

            for (int i = 0; i < mine; i++)
                in >> wantMine[i];

            for (int i = 0; i < terminal; i++)
                in >> wantTerminal[i];
        }

        in.close();
    }
    
    int print()
    {
        for (int i = 0; i < quantity; i++)
        {
            for (int j = 0; j < quantity; j++)
                std::cout << roads[i][j] << "\t";

            std::cout << std::endl;
        }

        for (int i = 0; i < mine; i++)
            std::cout << wantMine[i] << "\t";

        std::cout << std::endl;

        for (int i = 0; i < terminal; i++)
            std::cout << wantTerminal[i] << "\t";

        std::cout << std::endl;

        return 0;
    }
};