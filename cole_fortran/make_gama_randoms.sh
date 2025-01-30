#!/bin/bash

# Define the catalog identifiers
catalogs=("g09" "g12" "g15" "g23")

# Backup the original file
cp rancat_ascii_example.f90 rancat_ascii_example.bak

for cat in "${catalogs[@]}"; do
    echo "Updating the script for $cat"
    cp rancat_ascii_example.bak rancat_ascii_example.f90  # Restore original before modifying
    sed -i '' "s/catalogue.txt/${cat}_catalogue.txt/g" rancat_ascii_example.f90
    sed -i '' "s/rancat.txt/${cat}_rancat.txt/g" rancat_ascii_example.f90
    
    echo "Compiling and running Fortran"
    make
    ./rancat_ascii_example
    
    echo "Rewriting files for R"
    tail -n +7 "TestData/${cat}_rancat.txt" | awk '{print$2}' > "gama_${cat}_randoms.txt"
    sed -i '' 's/mag/z/g' "gama_${cat}_randoms.txt"
done

# Restore original file
echo "Cleaning up"
mv rancat_ascii_example.bak rancat_ascii_example.f90

echo "Randoms have been generated :)"
echo "Moving random catalogues to the Gama plotter folder..."
mv gama*randoms.txt ~/Desktop/GAMA_paper_plotter/