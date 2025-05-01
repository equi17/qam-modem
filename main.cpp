#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>
#include <vector>
#include <complex>
#include <fstream>

constexpr double SQRT_10 = 3.16227766017;
constexpr int BITS_PER_SYMBOL = 4;

constexpr std::complex<double> QAM16_SYMBOLS[16] = {
    {-3.0/SQRT_10,  3.0/SQRT_10}, // 0000
    {-3.0/SQRT_10,  1.0/SQRT_10}, // 0001
    {-3.0/SQRT_10, -1.0/SQRT_10}, // 0011
    {-3.0/SQRT_10, -3.0/SQRT_10}, // 0010

    {-1.0/SQRT_10,  3.0/SQRT_10}, // 0100
    {-1.0/SQRT_10,  1.0/SQRT_10}, // 0101
    {-1.0/SQRT_10, -1.0/SQRT_10}, // 0111
    {-1.0/SQRT_10, -3.0/SQRT_10}, // 0110

    { 1.0/SQRT_10,  3.0/SQRT_10}, // 1100
    { 1.0/SQRT_10,  1.0/SQRT_10}, // 1101
    { 1.0/SQRT_10, -1.0/SQRT_10}, // 1111
    { 1.0/SQRT_10, -3.0/SQRT_10}, // 1110

    { 3.0/SQRT_10,  3.0/SQRT_10}, // 1000
    { 3.0/SQRT_10,  1.0/SQRT_10}, // 1001
    { 3.0/SQRT_10, -1.0/SQRT_10}, // 1011
    { 3.0/SQRT_10, -3.0/SQRT_10}  // 1010
};

class QAM16Modulator {
public:
    std::vector<std::complex<double>> modulate(const std::vector<bool>& bits) {
        size_t padded_len = bits.size() + (BITS_PER_SYMBOL - (bits.size() % BITS_PER_SYMBOL)) % BITS_PER_SYMBOL;
        std::vector<std::complex<double>> symbols;
        symbols.reserve(padded_len / BITS_PER_SYMBOL);

        for (size_t i = 0; i < padded_len; i += BITS_PER_SYMBOL) {
            int idx = 0;
            for (size_t j = 0; j < BITS_PER_SYMBOL; ++j) {
                const bool bit = (i + j < bits.size()) ? bits[i + j] : false;
                idx |= (bit << (BITS_PER_SYMBOL - 1 - j));
            }
            symbols.push_back(QAM16_SYMBOLS[idx]);
        }

        return symbols;
    }
};

class QAM16Demodulator {
public:
    std::vector<bool> demodulate(const std::vector<std::complex<double>>& symbols) {
        std::vector<bool> bits;
        bits.reserve(symbols.size() * BITS_PER_SYMBOL);

        for (const auto& sym : symbols) {
            int closest_idx = find_closest_index(sym);
            for (int j = BITS_PER_SYMBOL - 1; j >= 0; --j) {
                bits.push_back((closest_idx >> j) & 1);
            }
        }

        return bits;
    }

private:
    int find_closest_index(const std::complex<double>& symbol) const {
        double min_dist = std::numeric_limits<double>::infinity();
        int closest_idx = 0;

        for (int i = 0; i < 16; ++i) {
            double dist = std::norm(symbol - QAM16_SYMBOLS[i]);
            if (dist < min_dist) {
                min_dist = dist;
                closest_idx = i;
            }
        }

        return closest_idx;
    }
};

class AWGN {
    std::mt19937 generator;
    std::normal_distribution<double> distribution;
    double current_stddev;

public:
    AWGN(double noise_stddev) 
        : generator(std::random_device{}()), 
          distribution(0.0, noise_stddev), 
          current_stddev(noise_stddev) {};

    void update_stddev(double new_stddev) {
        if (current_stddev != new_stddev) {
            distribution = std::normal_distribution<double>(0.0, new_stddev);
            current_stddev = new_stddev;
        }
    }

    void add_noise(std::vector<std::complex<double>>& symbols) {
        for (auto& sym : symbols) {
            sym += std::complex<double>(distribution(generator), distribution(generator));
        }
    }
};

class BitGenerator {
    std::mt19937 generator;
    std::bernoulli_distribution distribution;

public:
    BitGenerator() : generator(std::random_device{}()), distribution(0.5) {}

    std::vector<bool> generate(size_t num_bits) {
        std::vector<bool> bits;
        bits.reserve(num_bits);

        for (size_t i = 0; i < num_bits; ++i) {
            bits.push_back(distribution(generator));
        }

        return bits;
    }
};

int main() {
    size_t num_bits;
    double start_variance, end_variance, step_size;

    std::cout << "Введите количество бит: ";
    std::cin >> num_bits;

    std::cout << "Введите начальную дисперсию (σ²): ";
    std::cin >> start_variance;

    std::cout << "Введите конечную дисперсию (σ²): ";
    std::cin >> end_variance;

    std::cout << "Введите шаг по дисперсии (σ²): ";
    std::cin >> step_size;

    std::ofstream out_file("ber.csv");
    out_file << "noise_variance,ber\n";

    QAM16Modulator mod;
    QAM16Demodulator demod;
    AWGN noise(sqrt(start_variance));
    BitGenerator gen;

    std::vector<bool> bits = gen.generate(num_bits);
    std::vector<std::complex<double>> modulated_clean = mod.modulate(bits);

    for (double noise_variance = start_variance; noise_variance <= end_variance + 1e-9; noise_variance += step_size) {
        noise.update_stddev(sqrt(noise_variance));
        auto modulated = modulated_clean;
        noise.add_noise(modulated);
        std::vector<bool> out = demod.demodulate(modulated);

        int errors = 0;
        for (size_t j = 0; j < num_bits; ++j) {
            if (bits[j] != out[j]) ++errors;
        }

        double ber = static_cast<double>(errors) / num_bits;
        out_file << noise_variance << ',' << ber << '\n';
    }

    out_file.close();
    std::cout << "Данные записаны в ber.csv" << '\n';

    return 0;
}