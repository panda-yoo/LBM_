#include <iostream>
#include <fstream>
#include <array>
#include <filesystem>
// TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
#include "include/D1Q3.hpp"
#include "include/D2Q9.hpp"
// #include "include/D2Q9_cylinder.hpp"


int main() {
    // std::filesystem::create_directory("../data_LBM");
    d1q3_sim();
    // plane2d cyl;
    // cylinder_create(cyl);
    // show(cyl);
    // d2q9_sim();

    // sim_Cylinder_d2q9();
    return 0;
    // TIP See CLion help at <a href="https://www.jetbrains.com/help/clion/">jetbrains.com/help/clion/</a>. Also, you can try interactive lessons for CLion by selecting 'Help | Learn IDE Features' from the main menu.
}
