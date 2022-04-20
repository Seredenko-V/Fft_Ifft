#include <iostream>
#include <vector>
#include <iomanip> // округление
#include <complex>
#include <cmath>

#define PI 3.14159265358979323846
#define DIGITS_AFTER_DECIMAL_POINT 4

using namespace std;

vector <double> modulation(vector <uint32_t> bits, double Fn, double Fd) /* модуляция сигнала */
{
    double td = 1 / Fd; // шаг дискретизации по времени, сек
    uint32_t sample = uint32_t(Fd / Fn); // количество отсчетов на информационный символ
    vector <double> mod_code(sample * bits.size());
    for (size_t i = 0; i < bits.size(); i++)
    {
        //cout << endl << "Модуляция " << i << "-го бита: " << bits[i] << endl;
        for (size_t j = 0; j < sample; j++)
        {
            mod_code[j + sample * i] = sin(2 * PI * Fn * td * j + PI * bits[i]);
            //mod_code[j + sample * i] = sin(2 * PI * 1000 * td * j) + 0.5 * sin(2 * PI * 2000 * td * j + 3 * PI / 4);
            // + PI * bits[i]);
         //cout << fixed << showpoint << setprecision(DIGITS_AFTER_DECIMAL_POINT) << mod_code[j + sample * i] << endl;
        }
    }
    return mod_code;
}

vector<uint32_t> DecToBin(uint32_t dec_value, const uint32_t& size) {
    /* Перевод из 10-й в 2-ю СС. Взовращает массив бит */
    vector <uint32_t> bin_value;
    if (dec_value == 1)
        bin_value.push_back(1);
    if (dec_value == 0)
        bin_value.push_back(0);
    while (dec_value > 1) /* перевод в двоичную СС, результат записывается "задом наперед" */
    {
        bin_value.push_back(dec_value % 2);
        if (dec_value % 2 == 1)
            dec_value = dec_value - 1;
        dec_value = dec_value / 2;
        if (dec_value == 1)
            bin_value.push_back(dec_value % 2);
    }
    while (bin_value.size() < size)
        bin_value.push_back(0);
    vector<uint32_t> arr(bin_value.size()); // массив для корректного отображения числа в двоичной системе счисления
    for (int i = 0, j = bin_value.size() - 1; i < bin_value.size(), j >= 0; i++, j--)
    {
        arr[i] = bin_value[j];
    }
    for (int i = 0; i < bin_value.size(); i++)
    {
        bin_value[i] = arr[i];
    }
    return arr;
}

uint32_t BinToDec(const vector<uint32_t>& arr_bin, const uint32_t& size)
/* Перевод из 2-й в 10-ю СС. Взовращает число в 10-й СС */
{
    uint32_t dec_value = 0;
    for (size_t i = 0; i < size; i++)
    {
        dec_value += uint32_t(arr_bin[i] * pow(2, (size - 1 - i)));
    }
    return dec_value;
}

void inversion_arr(uint32_t* arr, uint32_t size)
/* Инверсия элементов массива. Инвертирует принимаемый массив */
{
    uint32_t buf;
    for (size_t i = 0, j = size - 1; i < j; i++, j--)
    {
        buf = arr[i];
        arr[i] = arr[j];
        arr[j] = buf;
    }
}

void CheckSize(vector<double>& sample) {
    /* Проверка соответствия размерности вектора степени двойки */
    if (log2(sample.size()) != int(log2(sample.size()))) {
        uint32_t new_size = uint32_t(pow(2, uint32_t(log2(sample.size())) + 1));
        sample.resize(new_size, 0); // дописываем нули к концу вектора
    }
}

vector<complex<double>> Fft(vector <double> sample) {
    /* Быстрое преобразование Фурье по основанию 2 */
    vector<complex<double>> sample_frequency(sample.size());
    CheckSize(sample);
    uint32_t size_bit = uint32_t(trunc(log2(sample.size() - 1) / log2(2)) + 1); // количество бит, занимаемых максимальным числом в массиве
    vector<vector<uint32_t>> bit_arr(sample.size(), vector<uint32_t>(size_bit)); // массив, хранящий двоичное представление порядковых номеров отсчетов
    for (size_t i = 0; i < sample.size(); ++i)
    {
        bit_arr[i] = DecToBin(i, size_bit);
    }
    vector<uint32_t> grouped_sample(sample.size()); // вектор, хранящий отсчеты в нужном порядке для БПФ
    for (size_t i = 0; i < sample.size(); ++i)
    {
        reverse(bit_arr[i].begin(), bit_arr[i].end());
        //inversion_arr(bit_arr[i], size_bit);
        grouped_sample[i] = BinToDec(bit_arr[i], size_bit);
    }
    // работа с комплексными числами:
    for (size_t i = 0; i < sample_frequency.size(); ++i) /* перевод вещественных остчетов в
                                                               комплексную форму согласно нужному порядку */
    {
        sample_frequency[i] = sample[grouped_sample[i]];
    }
    uint32_t number_of_tree_levels = uint32_t(log2(sample.size()) / log2(2)); /* количество каскадов дерева БПФ
                                                                  (сколько раз можно разделить пополам) */
    for (size_t i = 1; i <= number_of_tree_levels; ++i) // восходящие каскады БПФ
    {
        uint32_t m = uint32_t(pow(2, i));
        complex <double> Wm(cos(2 * PI / m), -sin(2 * PI / m)); /* поворачивающий множитель */
        for (size_t k = 0; k < sample.size(); k += m) // уменьшение числа ветвей
        {
            complex <double> W = 1;
            for (size_t j = 0; j < m / 2; ++j) // количество "перекрестий"
            {
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

vector<double> Ifft(const vector<complex<double>>& sample_frequency) {
    /* Обратное быстрое преобразование Фурье по основанию 2 */
    vector<double> time_sample(sample_frequency.size());
    uint32_t size_bit = uint32_t(trunc(log2(sample_frequency.size() - 1) / log2(2)) + 1); // количество бит, занимаемых максимальным числом в массиве
    vector<vector<uint32_t>> bit_arr(sample_frequency.size(), vector<uint32_t>(size_bit)); // массив, хранящий двоичное представление порядковых номеров отсчетов
    for (size_t i = 0; i < sample_frequency.size(); ++i)
    {
        bit_arr[i] = DecToBin(i, size_bit);
    }
    vector<uint32_t> grouped_sample(sample_frequency.size()); // вектор, хранящий отсчеты в нужном порядке для ОБПФ
    for (size_t i = 0; i < sample_frequency.size(); i++)
    {
        reverse(bit_arr[i].begin(), bit_arr[i].end());
        //inversion_arr(bit_arr[i], size_bit);
        grouped_sample[i] = BinToDec(bit_arr[i], size_bit);
    }

    // работа с комплексными числами:
    vector<complex<double>> grouped_sample_complex(sample_frequency.size());
    for (size_t i = 0; i < grouped_sample_complex.size(); ++i) /* перестановка комплексных остчетов
                                                               согласно нужному порядку для ОБПФ */
    {
        grouped_sample_complex[i] = sample_frequency[grouped_sample[i]];
    }

    uint32_t number_of_tree_levels = uint32_t(log2(sample_frequency.size()) / log2(2)); /* количество каскадов дерева ОБПФ
                                                                  (сколько раз можно разделить пополам) */
    for (size_t i = 1; i <= number_of_tree_levels; ++i) // восходящие каскады ОБПФ
    {
        uint32_t m = uint32_t(pow(2, i));
        complex <double> Wm(cos(2 * PI / m), sin(2 * PI / m)); /* поворачивающий множитель */
        for (size_t k = 0; k < sample_frequency.size(); k += m) // уменьшение числа ветвей
        {
            complex <double> W = 1;
            for (size_t j = 0; j < m / 2; ++j) // количество "перекрестий"
            {
                complex <double> odd = W * grouped_sample_complex[k + j + m / 2]; // нечетные
                complex <double> even = grouped_sample_complex[k + j]; // четные
                grouped_sample_complex[k + j] = even + odd;
                grouped_sample_complex[k + j + m / 2] = even - odd;
                W *= Wm;
            }
        }
    }
    for (size_t i = 0; i < sample_frequency.size(); ++i)
    {
        time_sample[i] = grouped_sample_complex[i].real() / grouped_sample_complex.size();
    }
    return time_sample;
}

int main()
{
    setlocale(LC_ALL, "rus");
    vector <uint32_t> inf_signal = { 0 };
    double Fn = 1000; // несущая частота, Гц
    double Fd = 8000; // частота дискретизации, Гц
    vector <double> mod_signal = modulation(inf_signal, Fn, Fd);
    for (size_t i = 0; i < mod_signal.size(); i++)
    {
        cout << fixed << showpoint << setprecision(DIGITS_AFTER_DECIMAL_POINT) << mod_signal[i] << endl;
    }
    cout << "========================" << endl;

    vector<complex<double>> result_fft = Fft(mod_signal);
    cout << "FFT" << endl;
    for (size_t i = 0; i < result_fft.size(); i++)
    {
        cout << result_fft[i] << endl;
    }
    cout << "========================" << endl;

    vector<double> time_sample = Ifft(result_fft);
    cout << "IFFT" << endl;
    for (size_t i = 0; i < time_sample.size(); i++)
    {
        cout << time_sample[i] << endl;
    }
    cout << "========================" << endl;

    return 0;
}