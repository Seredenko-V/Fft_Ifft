#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <iomanip>
#include <complex>
#include <cmath>
#include <algorithm>
#include <cassert>

const uint8_t DIGITS_AFTER_DECIMAL_POINT = 4u;

using namespace std;

// Абсолютная двухпозиционная фазовая модуляция. Возвращает вектор отсчетов сигнала
vector<double> AbsoluteTwoPositionPhaseModulation(const vector<uint32_t>& bits, const double Fn, const double Fd)  {
    double td = 1 / Fd; // шаг дискретизации по времени, сек
    uint32_t sample = uint32_t(Fd / Fn); // количество отсчетов на информационный символ
    vector<double> mod_code(sample * bits.size());
    for (size_t i = 0; i < bits.size(); i++) {
        for (size_t j = 0; j < sample; j++) {
            mod_code[j + sample * i] = sin(2 * M_PI * Fn * td * j + M_PI * bits[i]);
        }
    }
    return mod_code;
}

namespace convert_digits {
    // Перевод из 10-й в 2-ю СС. Взовращает вектор бит
    vector<uint32_t> DecToBin(uint32_t dec_value, const uint32_t size) {
        vector<uint32_t> bin_value;
        if (dec_value == 1) {
            bin_value.push_back(1);
        } else if (dec_value == 0) {
            bin_value.push_back(0);
        }

        // перевод в двоичную СС, результат записывается "задом наперед"
        while (dec_value > 1) {
            bin_value.push_back(dec_value % 2);
            if (dec_value % 2 == 1) {
                dec_value = dec_value - 1;
            }
            dec_value = dec_value / 2;
            if (dec_value == 1) {
                bin_value.push_back(dec_value % 2);
            }
        }

        bin_value.resize(size, 0); // дописываем нули к концу вектора
        reverse(bin_value.begin(), bin_value.end());
        return bin_value;
    }

    // Перевод из 2-й в 10-ю СС. Взовращает число в 10-й СС
    uint32_t BinToDec(const vector<uint32_t>& arr_bin, const uint32_t size) {
        uint32_t dec_value = 0;
        for (size_t i = 0; i < size; i++) {
            dec_value += uint32_t(arr_bin[i] * pow(2, (size - 1 - i)));
        }
        return dec_value;
    }
} // namespace convert_digits

namespace fft_ifft {

    // Проверка соответствия размерности вектора степени двойки
    void CheckSize(vector<double>& sample) {
        if (log2(sample.size()) != int(log2(sample.size()))) {
            uint32_t new_size = uint32_t(pow(2, uint32_t(log2(sample.size())) + 1));
            sample.resize(new_size, 0); // дописываем нули к концу вектора
        }
    }

    // Быстрое преобразование Фурье по основанию 2
    vector<complex<double>> Fft(vector<double> sample) {
        using namespace convert_digits;
        vector<complex<double>> sample_frequency(sample.size());
        CheckSize(sample);

        uint32_t size_bit = uint32_t(trunc(log2(sample.size() - 1) / log2(2)) + 1); // количество бит, занимаемых максимальным числом в массиве
        vector<vector<uint32_t>> bit_arr(sample.size(), vector<uint32_t>(size_bit)); // массив, хранящий двоичное представление порядковых номеров отсчетов
        for (size_t i = 0; i < sample.size(); ++i) {
            bit_arr[i] = DecToBin(static_cast<uint32_t>(i), size_bit);
        }

        vector<uint32_t> grouped_sample(sample.size()); // вектор, хранящий отсчеты в нужном порядке для БПФ
        for (size_t i = 0; i < sample.size(); ++i) {
            reverse(bit_arr[i].begin(), bit_arr[i].end());
            grouped_sample[i] = BinToDec(bit_arr[i], size_bit);
        }

        // работа с комплексными числами:
        // перевод вещественных остчетов в комплексную форму согласно нужному порядку
        for (size_t i = 0; i < sample_frequency.size(); ++i) {
            sample_frequency[i] = sample[grouped_sample[i]];
        }
        uint32_t number_of_tree_levels = uint32_t(log2(sample.size()) / log2(2)); /* количество каскадов дерева БПФ
                                                                      (сколько раз можно разделить пополам) */
                                                                      // восходящие каскады БПФ
        for (size_t i = 1; i <= number_of_tree_levels; ++i) {
            uint32_t m = uint32_t(pow(2, i));
            complex <double> Wm(cos(2 * M_PI / m), -sin(2 * M_PI / m)); /* поворачивающий множитель */
            // уменьшение числа ветвей
            for (size_t k = 0; k < sample.size(); k += m) {
                complex <double> W = 1;
                // количество "перекрестий"
                for (size_t j = 0; j < m / 2; ++j) {
                    complex <double> odd = W * sample_frequency[k + j + m / 2]; // нечетные
                    complex <double> even = sample_frequency[k + j]; // четные
                    sample_frequency[k + j] = even + odd;
                    sample_frequency[k + j + m / 2] = even - odd;
                    W *= Wm;
                }
            }
        }
        return sample_frequency;
    }

    // Обратное быстрое преобразование Фурье по основанию 2
    vector<double> Ifft(const vector<complex<double>>& sample_frequency) {
        using namespace convert_digits;
        vector<double> time_sample(sample_frequency.size());

        // количество бит, занимаемых максимальным числом в массиве
        uint32_t size_bit = uint32_t(trunc(log2(sample_frequency.size() - 1) / log2(2)) + 1);

        // массив, хранящий двоичное представление порядковых номеров отсчетов
        vector<vector<uint32_t>> bit_arr(sample_frequency.size(), vector<uint32_t>(size_bit));
        for (size_t i = 0; i < sample_frequency.size(); ++i) {
            bit_arr[i] = DecToBin(static_cast<uint32_t>(i), size_bit);
        }

        // вектор, хранящий отсчеты в нужном порядке для ОБПФ
        vector<uint32_t> grouped_sample(sample_frequency.size());
        for (size_t i = 0; i < sample_frequency.size(); i++) {
            reverse(bit_arr[i].begin(), bit_arr[i].end());
            grouped_sample[i] = BinToDec(bit_arr[i], size_bit);
        }

        // работа с комплексными числами:
        vector<complex<double>> grouped_sample_complex(sample_frequency.size());

        // перестановка комплексных остчетов согласно нужному порядку для ОБПФ
        for (size_t i = 0; i < grouped_sample_complex.size(); ++i) {
            grouped_sample_complex[i] = sample_frequency[grouped_sample[i]];
        }

        // количество каскадов дерева ОБПФ (сколько раз можно разделить пополам)
        uint32_t number_of_tree_levels = uint32_t(log2(sample_frequency.size()) / log2(2));

        // восходящие каскады ОБПФ
        for (size_t i = 1; i <= number_of_tree_levels; ++i) {
            uint32_t m = uint32_t(pow(2, i));
            complex <double> Wm(cos(2 * M_PI / m), sin(2 * M_PI / m)); // поворачивающий множитель
            // уменьшение числа ветвей
            for (size_t k = 0; k < sample_frequency.size(); k += m) {
                complex <double> W = 1;
                // количество "перекрестий"
                for (size_t j = 0; j < m / 2; ++j) {
                    complex <double> odd = W * grouped_sample_complex[k + j + m / 2]; // нечетные
                    complex <double> even = grouped_sample_complex[k + j]; // четные
                    grouped_sample_complex[k + j] = even + odd;
                    grouped_sample_complex[k + j + m / 2] = even - odd;
                    W *= Wm;
                }
            }
        }
        for (size_t i = 0; i < sample_frequency.size(); ++i) {
            time_sample[i] = grouped_sample_complex[i].real() / grouped_sample_complex.size();
        }
        return time_sample;
    }
} // namespace FftIfft

template <typename Type>
void PrintVector(ostream& out, const vector<Type>& print_vector) {
    out << fixed << showpoint << setprecision(DIGITS_AFTER_DECIMAL_POINT);
    for (const Type& element : print_vector) {
        out << element << endl;
    }
    out << "========================" << endl;
}

void Test(bool do_print) {
    using namespace fft_ifft;
    const double DELTA = 1e-6;
    vector<uint32_t> inf_signal = { 0 };
    double Fn = 1000; // несущая частота, Гц
    double Fd = 8000; // частота дискретизации, Гц
    vector<double> mod_signal = AbsoluteTwoPositionPhaseModulation(inf_signal, Fn, Fd);
    vector<complex<double>> result_fft = Fft(mod_signal);
    vector<double> time_sample = Ifft(result_fft);
    assert(mod_signal.size() == time_sample.size());
    for (size_t i = 0; i < mod_signal.size(); ++i) {
        assert(abs(mod_signal[i] - time_sample[i]) <= DELTA);
    }

    if (do_print) {
        cout << "Etalon signal:" << endl;
        PrintVector(cout, mod_signal);
        cout << "Result Fft:" << endl;
        PrintVector(cout, result_fft);
        cout << "Result Ifft:" << endl;
        PrintVector(cout, time_sample);
    }
}

int main() {
    bool do_print_result = true;
    Test(do_print_result);
    return 0;
}