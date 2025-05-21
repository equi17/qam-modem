#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>
#include <vector>
#include <complex>
#include <fstream>

constexpr double SQRT_2 = 1.41421356237;
constexpr double SQRT_10 = 3.16227766017;
constexpr double SQRT_42 = 6.48074069841;

constexpr double AVG_SYMBOL_ENERGY = 1.0;

constexpr std::complex<double> QPSK_SYMBOLS[4] = {
    { 1.0/SQRT_2,  1.0/SQRT_2}, // 00
    { 1.0/SQRT_2, -1.0/SQRT_2}, // 01
    {-1.0/SQRT_2,  1.0/SQRT_2}, // 11
    {-1.0/SQRT_2, -1.0/SQRT_2}  // 10
};

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

constexpr std::complex<double> QAM64_SYMBOLS[64] = {
    {-7.0/SQRT_42,  7.0/SQRT_42}, {-7.0/SQRT_42,  5.0/SQRT_42}, 
    {-7.0/SQRT_42,  3.0/SQRT_42}, {-7.0/SQRT_42,  1.0/SQRT_42},
    {-7.0/SQRT_42, -1.0/SQRT_42}, {-7.0/SQRT_42, -3.0/SQRT_42},
    {-7.0/SQRT_42, -5.0/SQRT_42}, {-7.0/SQRT_42, -7.0/SQRT_42},

    {-5.0/SQRT_42,  7.0/SQRT_42}, {-5.0/SQRT_42,  5.0/SQRT_42},
    {-5.0/SQRT_42,  3.0/SQRT_42}, {-5.0/SQRT_42,  1.0/SQRT_42},
    {-5.0/SQRT_42, -1.0/SQRT_42}, {-5.0/SQRT_42, -3.0/SQRT_42},
    {-5.0/SQRT_42, -5.0/SQRT_42}, {-5.0/SQRT_42, -7.0/SQRT_42},

    {-3.0/SQRT_42,  7.0/SQRT_42}, {-3.0/SQRT_42,  5.0/SQRT_42},
    {-3.0/SQRT_42,  3.0/SQRT_42}, {-3.0/SQRT_42,  1.0/SQRT_42},
    {-3.0/SQRT_42, -1.0/SQRT_42}, {-3.0/SQRT_42, -3.0/SQRT_42},
    {-3.0/SQRT_42, -5.0/SQRT_42}, {-3.0/SQRT_42, -7.0/SQRT_42},

    {-1.0/SQRT_42,  7.0/SQRT_42}, {-1.0/SQRT_42,  5.0/SQRT_42},
    {-1.0/SQRT_42,  3.0/SQRT_42}, {-1.0/SQRT_42,  1.0/SQRT_42},
    {-1.0/SQRT_42, -1.0/SQRT_42}, {-1.0/SQRT_42, -3.0/SQRT_42},
    {-1.0/SQRT_42, -5.0/SQRT_42}, {-1.0/SQRT_42, -7.0/SQRT_42},

    { 1.0/SQRT_42,  7.0/SQRT_42}, { 1.0/SQRT_42,  5.0/SQRT_42},
    { 1.0/SQRT_42,  3.0/SQRT_42}, { 1.0/SQRT_42,  1.0/SQRT_42},
    { 1.0/SQRT_42, -1.0/SQRT_42}, { 1.0/SQRT_42, -3.0/SQRT_42},
    { 1.0/SQRT_42, -5.0/SQRT_42}, { 1.0/SQRT_42, -7.0/SQRT_42},

    { 3.0/SQRT_42,  7.0/SQRT_42}, { 3.0/SQRT_42,  5.0/SQRT_42},
    { 3.0/SQRT_42,  3.0/SQRT_42}, { 3.0/SQRT_42,  1.0/SQRT_42},
    { 3.0/SQRT_42, -1.0/SQRT_42}, { 3.0/SQRT_42, -3.0/SQRT_42},
    { 3.0/SQRT_42, -5.0/SQRT_42}, { 3.0/SQRT_42, -7.0/SQRT_42},

    { 5.0/SQRT_42,  7.0/SQRT_42}, { 5.0/SQRT_42,  5.0/SQRT_42},
    { 5.0/SQRT_42,  3.0/SQRT_42}, { 5.0/SQRT_42,  1.0/SQRT_42},
    { 5.0/SQRT_42, -1.0/SQRT_42}, { 5.0/SQRT_42, -3.0/SQRT_42},
    { 5.0/SQRT_42, -5.0/SQRT_42}, { 5.0/SQRT_42, -7.0/SQRT_42},

    { 7.0/SQRT_42,  7.0/SQRT_42}, { 7.0/SQRT_42,  5.0/SQRT_42},
    { 7.0/SQRT_42,  3.0/SQRT_42}, { 7.0/SQRT_42,  1.0/SQRT_42},
    { 7.0/SQRT_42, -1.0/SQRT_42}, { 7.0/SQRT_42, -3.0/SQRT_42},
    { 7.0/SQRT_42, -5.0/SQRT_42}, { 7.0/SQRT_42, -7.0/SQRT_42}
};

enum class ModulationType { QPSK, QAM16, QAM64 };

class QAMModulator {
    ModulationType type;
    int bits_per_symbol;
    std::vector<std::complex<double>> constellation;

public:
    QAMModulator(ModulationType type) : type(type) {
        switch (type) {
            case ModulationType::QPSK:
                bits_per_symbol = 2;
                constellation.reserve(4);
                std::copy(&QPSK_SYMBOLS[0], &QPSK_SYMBOLS[4], std::back_inserter(constellation));
                break;
            case ModulationType::QAM16:
                bits_per_symbol = 4;
                constellation.reserve(16);
                std::copy(&QAM16_SYMBOLS[0], &QAM16_SYMBOLS[16], std::back_inserter(constellation));
                break;
            case ModulationType::QAM64:
                bits_per_symbol = 6;
                constellation.reserve(64);
                std::copy(&QAM64_SYMBOLS[0], &QAM64_SYMBOLS[64], std::back_inserter(constellation));
                break;
        }
    }

    std::vector<std::complex<double>> modulate(const std::vector<bool>& bits) {
        size_t padded_len = bits.size() + (bits_per_symbol - (bits.size() % bits_per_symbol)) % bits_per_symbol;
        std::vector<std::complex<double>> symbols;
        symbols.reserve(padded_len / bits_per_symbol);

        for (size_t i = 0; i < padded_len; i += bits_per_symbol) {
            int idx = 0;
            for (size_t j = 0; j < bits_per_symbol; ++j) {
                const bool bit = (i + j < bits.size()) ? bits[i + j] : false;
                idx |= (bit << (bits_per_symbol - 1 - j));
            }
            symbols.push_back(constellation[idx]);
        }

        return symbols;
    }
};

class QAMDemodulator {
ModulationType type;
int bits_per_symbol;
std::vector<std::complex<double>> constellation;

public:
    QAMDemodulator(ModulationType type) : type(type) {
        switch (type) {
            case ModulationType::QPSK:
                bits_per_symbol = 2;
                constellation.reserve(4);
                std::copy(&QPSK_SYMBOLS[0], &QPSK_SYMBOLS[4], std::back_inserter(constellation));
                break;
            case ModulationType::QAM16:
                bits_per_symbol = 4;
                constellation.reserve(16);
                std::copy(&QAM16_SYMBOLS[0], &QAM16_SYMBOLS[16], std::back_inserter(constellation));
                break;
            case ModulationType::QAM64:
                bits_per_symbol = 6;
                constellation.reserve(64);
                std::copy(&QAM64_SYMBOLS[0], &QAM64_SYMBOLS[64], std::back_inserter(constellation));
                break;
        }
    }
    
    std::vector<bool> demodulate(const std::vector<std::complex<double>>& symbols) {
        std::vector<bool> bits;
        bits.reserve(symbols.size() * bits_per_symbol);

        for (const auto& sym : symbols) {
            int closest_idx = find_closest_index(sym);
            for (int j = bits_per_symbol - 1; j >= 0; --j) {
                bits.push_back((closest_idx >> j) & 1);
            }
        }

        return bits;
    }

private:
    int find_closest_index(const std::complex<double>& symbol) const {
        double min_dist = std::numeric_limits<double>::infinity();
        int closest_idx = 0;

        for (int i = 0; i < constellation.size(); ++i) {
            double dist = std::norm(symbol - constellation[i]);
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
    double start_snr_db, end_snr_db, snr_step_db;
    int modulation_type;

    std::cout << "Выберите тип модуляции (0-QPSK, 1-QAM16, 2-QAM64): ";
    std::cin >> modulation_type;
    
    ModulationType mod_type;
    switch (modulation_type) {
        case 0: mod_type = ModulationType::QPSK; break;
        case 1: mod_type = ModulationType::QAM16; break;
        case 2: mod_type = ModulationType::QAM64; break;
        default:
            std::cerr << "Неверный тип модуляции\n";
            return 1;
    }

    std::cout << "Введите количество бит: ";
    std::cin >> num_bits;

    std::cout << "Введите начальное SNR (dB): ";
    std::cin >> start_snr_db;
    
    std::cout << "Введите конечное SNR (dB): ";
    std::cin >> end_snr_db;
    
    std::cout << "Введите шаг по SNR (dB): ";
    std::cin >> snr_step_db;

    std::string filename;
    switch (mod_type) {
        case ModulationType::QPSK: filename = "qpsk_ber.csv"; break;
        case ModulationType::QAM16: filename = "qam16_ber.csv"; break;
        case ModulationType::QAM64: filename = "qam64_ber.csv"; break;
    }

    double Es = AVG_SYMBOL_ENERGY;

    std::ofstream out_file(filename);
    out_file << "snr_db,ber\n";

    QAMModulator mod(mod_type);
    QAMDemodulator demod(mod_type);
    AWGN noise(0.1);
    BitGenerator gen;

    std::vector<bool> bits = gen.generate(num_bits);
    std::vector<std::complex<double>> modulated_clean = mod.modulate(bits);

    for (double snr_db = start_snr_db; snr_db <= end_snr_db + 1e-9; snr_db += snr_step_db) {
        double snr_linear = pow(10.0, snr_db / 10.0);
        double noise_variance = Es / (2.0 * snr_linear);
        
        noise.update_stddev(sqrt(noise_variance));
        auto modulated = modulated_clean;
        noise.add_noise(modulated);
        
        std::vector<bool> out = demod.demodulate(modulated);

        int errors = 0;
        for (size_t j = 0; j < bits.size(); ++j) {
            if (bits[j] != out[j]) ++errors;
        }

        double ber = static_cast<double>(errors) / bits.size();
        out_file << snr_db << ',' << ber << '\n';
    }

    out_file.close();
    std::cout << "Данные записаны в " << filename << '\n';

    return 0;
}