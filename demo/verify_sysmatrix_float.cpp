#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <cstdint>
#include <cstring> // Required for std::memcmp

// 颜色输出宏
#define RED "\033[31m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define RESET "\033[0m"

struct FileStats {
    long long fileSize;
    long long dataSegmentSize; // 除去尾部浮点数的大小
};

FileStats getFileStats(const std::string& filename, int Nx, int Ny) {
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    if (!file) {
        std::cerr << RED << "Error: Cannot open file " << filename << RESET << std::endl;
        exit(1);
    }
    long long size = file.tellg();
    long long floatSize = (long long)Nx * Ny * sizeof(float);
    
    if (size < floatSize) {
        std::cerr << RED << "Error: File size (" << size << ") is smaller than expected float map size (" << floatSize << ")" << RESET << std::endl;
        exit(1);
    }
    
    return {size, size - floatSize};
}

void compare_binary_body(const std::string& f1_name, const std::string& f2_name, long long limit) {
    std::cout << "--- Phase 1: Checking Binary Consistency (Structure & Quantized Data) ---" << std::endl;
    
    std::ifstream f1(f1_name, std::ios::binary);
    std::ifstream f2(f2_name, std::ios::binary);
    
    const size_t buffer_size = 1024 * 1024 * 16; // 16MB buffer
    std::vector<char> buf1(buffer_size);
    std::vector<char> buf2(buffer_size);
    
    long long processed = 0;
    while (processed < limit) {
        long long to_read = std::min((long long)buffer_size, limit - processed);
        f1.read(buf1.data(), to_read);
        f2.read(buf2.data(), to_read);
        
        if (std::memcmp(buf1.data(), buf2.data(), to_read) != 0) {
            // Find exact error location
            for (long long i = 0; i < to_read; ++i) {
                if (buf1[i] != buf2[i]) {
                    std::cerr << RED << "[FAIL] Binary mismatch at byte offset: " << processed + i << RESET << std::endl;
                    std::cerr << "File 1: " << std::hex << (int)(uint8_t)buf1[i] 
                              << ", File 2: " << (int)(uint8_t)buf2[i] << std::dec << std::endl;
                    exit(1);
                }
            }
        }
        processed += to_read;
        if (processed % (1024*1024*100) == 0) {
            std::cout << "\rChecked " << processed / (1024*1024) << " MB..." << std::flush;
        }
    }
    std::cout << "\r" << GREEN << "[PASS] Binary body (" << limit << " bytes) is identical." << RESET << std::endl;
}

void compare_float_tail(const std::string& f1_name, const std::string& f2_name, int Nx, int Ny) {
    std::cout << "\n--- Phase 2: Analyzing Floating Point Differences (Aval_max_ptr) ---" << std::endl;
    
    long long float_count = (long long)Nx * Ny;
    long long offset_from_end = float_count * sizeof(float);
    
    std::ifstream f1(f1_name, std::ios::binary);
    std::ifstream f2(f2_name, std::ios::binary);
    
    // Seek to the float section (End of file - Nx*Ny*4)
    f1.seekg(-offset_from_end, std::ios::end);
    f2.seekg(-offset_from_end, std::ios::end);
    
    std::vector<float> data1(float_count);
    std::vector<float> data2(float_count);
    
    f1.read(reinterpret_cast<char*>(data1.data()), float_count * sizeof(float));
    f2.read(reinterpret_cast<char*>(data2.data()), float_count * sizeof(float));
    
    double max_diff = 0.0;
    double sum_sq_diff = 0.0;
    long long diff_count = 0;
    int print_limit = 10;
    
    std::cout << "Comparing " << float_count << " float values..." << std::endl;
    
    for (long long i = 0; i < float_count; ++i) {
        float v1 = data1[i];
        float v2 = data2[i];
        float diff = std::abs(v1 - v2);
        
        if (diff > 0.0) {
            if (diff > max_diff) max_diff = diff;
            sum_sq_diff += (double)diff * diff;
            diff_count++;
            
            if (print_limit > 0 && diff > 1e-5) { // Only print significant errors
                std::cout << YELLOW << "Diff at index " << i << " (" << i/Nx << "," << i%Nx << "): "
                          << v1 << " vs " << v2 << " | Delta: " << diff << RESET << std::endl;
                print_limit--;
            }
        }
    }
    
    double mse = sum_sq_diff / float_count;
    double rmse = std::sqrt(mse);
    
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "Total Floats:   " << float_count << std::endl;
    std::cout << "Differing Vals: " << diff_count << " (" << (double)diff_count/float_count*100.0 << "%)" << std::endl;
    std::cout << "Max Difference: " << (max_diff > 1e-4 ? RED : GREEN) << max_diff << RESET << std::endl;
    std::cout << "MSE:            " << mse << std::endl;
    std::cout << "RMSE:           " << rmse << std::endl;
    
    if (max_diff < 1e-5) {
        std::cout << GREEN << "[SUCCESS] Floating point values match within tolerance." << RESET << std::endl;
    } else {
        std::cout << YELLOW << "[WARNING] Floating point differences detected." << RESET << std::endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cout << "Usage: " << argv[0] << " <file1> <file2> <Nx> <Ny>" << std::endl;
        std::cout << "Example: " << argv[0] << " orig.2Dsvmatrix stream.2Dsvmatrix 512 512" << std::endl;
        return 1;
    }
    
    std::string f1 = argv[1];
    std::string f2 = argv[2];
    int Nx = std::atoi(argv[3]);
    int Ny = std::atoi(argv[4]);
    
    std::cout << "Comparing File: " << f1 << std::endl;
    std::cout << "          With: " << f2 << std::endl;
    std::cout << "Dimensions:     " << Nx << " x " << Ny << std::endl;
    
    FileStats s1 = getFileStats(f1, Nx, Ny);
    FileStats s2 = getFileStats(f2, Nx, Ny);
    
    if (s1.fileSize != s2.fileSize) {
        std::cerr << RED << "[FAIL] File sizes differ!" << RESET << std::endl;
        std::cerr << f1 << ": " << s1.fileSize << " bytes" << std::endl;
        std::cerr << f2 << ": " << s2.fileSize << " bytes" << std::endl;
        std::cerr << "Diff: " << std::abs(s1.fileSize - s2.fileSize) << " bytes" << std::endl;
        return 1;
    }
    
    // 1. 比较二进制主体 (除了最后的 Nx*Ny 个浮点数)
    compare_binary_body(f1, f2, s1.dataSegmentSize);
    
    // 2. 比较尾部的浮点数
    compare_float_tail(f1, f2, Nx, Ny);
    
    return 0;
}
