//
// Created by au642261 on 4/24/25.
//

#ifndef GRID2D_H
#define GRID2D_H

#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>

template<typename T>
class Grid2D {
public:
    Grid2D() = default;

    Grid2D(size_t nRows, size_t nCols)
        : nx(nRows), ny(nCols), data(nx * ny) {}

    // Returns reference (not const allowing assignment)
    T& operator()(size_t i, size_t j) {
        //check_bounds(i, j);
        return data[i * ny + j];
    }

    // Returns const reference (read-only)
    const T& operator()(size_t i, size_t j) const {
        //check_bounds(i, j);
        return data[i * ny + j];
    }

    // Same functions as above, but here with bounds checking
    T& at(size_t i, size_t j) {
        check_bounds(i,j);
        return data[i * ny + j];
    }
    const T& at(size_t i, size_t j) const {
        check_bounds(i,j);
        return data[i * ny + j];
    }

    // Method for filling the data vector with one value
    void fill(T value) {
        std::fill(data.begin(), data.end(), value);
    }

    // Method to transpose the grid using a temporary array
    void transpose_inplace() {
        // Create a temporary grid to hold the transposed data
        std::vector<T> temp(ny * nx); // Transposed dimensions (ny rows, nx columns)

        // Copy the data into the temporary array, transposing the rows and columns
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                temp[j * nx + i] = (*this)(i, j); // swap row and column indices
            }
        }

        // Now assign the transposed data back to the original grid
        data = std::move(temp);  // Move the data from temp to data

        // Swap the dimensions
        std::swap(nx, ny);
    }

    // Method to return pointer to raw data?
    T* raw_data() { return data.data(); }
    const T* raw_data() const { return data.data(); }

    // Method for returning reference to data?
    //std::vector<T>& dataReference() {return &data;}
    const std::vector<T>& dataReference() const {return data;}

    // Method for printing out the data grid
    void printGrid() const {
        for (size_t i=0; i<nx; i++) {
            for (size_t j=0; j<ny; j++) {
                std::cout << this->operator()(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }

    // Method to save to file
    void saveToFile(const std::string& filename) const {
        std::ofstream out(filename);
        if (!out) {
            std::cerr << "Error opening file!\n";
            return;
        }

        // Set the precision to 10 digits
        out << std::fixed << std::setprecision(14); // fixed-point notation with 10 digits

        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                out << this->operator()(i, j);
                if (j < ny - 1) out << " ";
            }
            out << "\n";
        }
    }

    size_t nRows() const { return nx; }
    size_t nCols() const { return ny; }


private:
    size_t nx, ny;
    std::vector<T> data;

    void check_bounds(size_t i, size_t j) const {
        if (i >= nx || j >= ny) {
            throw std::out_of_range("Grid2D index out of bounds");
        }
    }
};

#endif //GRID2D_H
