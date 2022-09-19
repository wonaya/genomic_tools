// my first program in C++

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

std::ifstream file("test.out");
std::string line;

while(std::getline(file, line))
{
    std::stringstream   linestream(line);
    std::string         data;
    int                 val1;
    int                 val2;

    std::getline(linestream, data, '\t');  // read up-to the first tab (discard tab).

    linestream >> val1 >> val2;
}