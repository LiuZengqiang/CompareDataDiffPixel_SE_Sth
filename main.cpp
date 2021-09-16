#include <stdio.h>
#include <iostream>
#include <tuple>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>
#include <assert.h>

using namespace std;
const double Epsilon_1 = 1e-6;
const double Epsilon_2 = 2e-3;

void process(float percent);

// R^2 = 1 - (SS_res)/(SS_tot)
float getCoefficientOfDetermination(vector<vector<float>>& measure, vector<vector<float>>& ground_truth);

void done(float origin_pixel_size, vector<int>& pixel_scale_factors, vector<float>& pixel_sizes, string input_file_path, string output_file_path, float heliostat_area, float receiver_height, float receiver_width);

vector<vector<float>> getDataMatrix(fstream& file);

vector<vector<float>> getSubSampledMatrix(vector<vector<float>>& origin, int dh = 1, int dw = 1);

void split(string str, vector<float>& v, string spacer);

void splitStr(string str, vector<string>& v, string spacer);

void
combineResults(vector<string>& input_file_lists, vector<string> titles, string target_title, string output_file_path);


int main(int arc, char** argv) {
    float origin_pixel = 0.01;
    vector<int> pixel_scale_factor;
    vector<float> pixel_sizes;

    for (int i = 2; i <= 20; i++) {
        pixel_scale_factor.push_back(i);
        pixel_sizes.push_back(origin_pixel * i);
    }

    float receiver_height = 20.0;
    float receiver_width = 20.0;
    float heliostat_area = 4.0f * 3.2f;

    vector<int > distance;
    for (int i = 200; i < 1200; i += 50) {
        distance.push_back(i);
    }
    string output_file_path_pre = "../Output/out_";
    string output_file_path_suf = ".csv";
    for (int dis : distance) {
        string input_file_path_pre = "/home/sth/CLionProjects/SolarEnergy_Chier/OutputFiles/GroundTruth_6282/Altitude_90/Normal_1/Dis_" + to_string(dis) + "/";
        string input_file_path_suf = "before_smooth_d_2048_r_128.csv";

        string output_file_path = output_file_path_pre + to_string(dis) + output_file_path_suf;
        done(origin_pixel, pixel_scale_factor, pixel_sizes, input_file_path_pre + input_file_path_suf, output_file_path, heliostat_area, receiver_height, receiver_width);
    }

    string combine_output_path_pre = "../Output/out_all_";
    string combine_output_path_suf = ".csv";

    string title = "dis_";
    // ----- 1. Calculate the compare result -----
    // 输入文件夹列表，每个文件夹输出一个结果(一个excel表格)
    // // ----- 2. Combine the compare results -----
    vector<string> target = { "rms", "max", "r2", "sum", "rms_per_area", "total_energy_err", "abs_energy_err" };
    vector<string> output_file_paths = { "../Output/out_200.csv", "../Output/out_500.csv","../Output/out_1000.csv" };
    for (int dis : distance) {
        output_file_paths.push_back(output_file_path_pre + to_string(dis) + output_file_path_suf);
    }

    string title_pre = title;
    string title_suf = "";
    vector<string> titles;
    for (int dis : distance) {
        titles.push_back(title_pre + to_string(dis) + title_suf);
    }
    for (string t : target) {
        string combine_output_path = combine_output_path_pre + t + combine_output_path_suf;
        combineResults(output_file_paths, titles, t, combine_output_path);
    }
    return 0;
}


void process(float percent) {
    percent = min(1.0f, percent);
    percent = max(0.0f, percent);

    int num = static_cast<int>(percent / 0.1f);
    cout << "[";
    for (int i = 0; i < num; i++) {
        cout << "=";
    }
    cout << ">";
    for (int i = num; i < 10; i++) {
        cout << " ";
    }

    cout << static_cast<int>(percent * 100) << "%" << endl;

}

/**
 * Coefficient of determination is equal to
 * @param measure
 * @param ground_truth
 * @return
 */
float getCoefficientOfDetermination(vector<vector<float>>& measure, vector<vector<float>>& ground_truth) {
    float SS_res = 0.0f;
    float SS_tot = 0.0f;
    float average_y = 0.0f;
    int ground_truth_row = ground_truth.size();
    int ground_truth_col = ground_truth[0].size();
    int measure_row = measure.size();
    int measure_col = measure[0].size();
    int factor_row = ground_truth_row / measure.size();
    int factor_col = ground_truth_col / measure[0].size();

    for (int row = 0; row < ground_truth_row; row++) {
        for (int col = 0; col < ground_truth_col; col++) {
            int r = min(measure_row - 1, row / factor_row);
            int c = min(measure_col - 1, col / factor_col);
            SS_res += pow(measure[r][c] - ground_truth[row][col], 2);
            average_y += measure[r][c];
        }
    }

    average_y /= ground_truth_row * ground_truth_col;
    for (int row = 0; row < ground_truth_row; row++) {
        for (int col = 0; col < ground_truth_col; col++) {
            int r = min(measure_row - 1, row / factor_row);
            int c = min(measure_col - 1, col / factor_col);
            SS_tot += pow(measure[r][c] - average_y, 2);
        }
    }
    return 1.0f - SS_res / SS_tot;
}
// void done(float origin_pixel_size, vector<int>& pixel_scale_factors, vector<float>& pixel_sizes, string input_file_path, string output_file_path, float heliostat_area, float receiver_height, float receiver_width);

void done(float origin_pixel_size, vector<int>& pixel_scale_factors, vector<float>& pixel_sizes, string input_file_path, string output_file_path, float heliostat_area, float receiver_height, float receiver_width) {

    fstream output_file(output_file_path, fstream::out);

    output_file
        << "pixel_scale_factor,pixel_size,max,sum,rms,ave,pixel,r2,rms_per_area,total_energy_err,abs_energy_err"
        << endl;
    cout << "Ground truth file path:" << input_file_path << endl;

    fstream ground_truth_file(input_file_path.c_str());
    if (!ground_truth_file.good()) {
        cerr << "Ground truth file " << input_file_path << " opened fail." << endl;
        return;
    }
    vector<vector<float>> ground_truth = getDataMatrix(ground_truth_file);

    // get result data size (h,w)
    int ground_truth_row = ground_truth.size();
    int ground_truth_column = ground_truth[0].size();
    float ground_truth_energy = 0.0f;

    float receiver_pixel_height = receiver_height / ground_truth_row;
    float receiver_pixel_width = receiver_width / ground_truth_column;

    float receiver_pixel_area = receiver_pixel_height * receiver_pixel_width;

    float ground_truth_average = 0.0f;


    for (int row = 0; row < ground_truth_row; row++) {
        for (int col = 0; col < ground_truth_column; col++) {
            ground_truth_energy += (receiver_pixel_area * ground_truth[row][col]);
            ground_truth_average += ground_truth[row][col];
        }
    }

    ground_truth_average /= (ground_truth_row * ground_truth_column);


    // abs_energy
    float max_dif, sum_dif, rms_dif, ave_dif, R2, rms_per_area, total_energy_err, abs_energy_err;
    int pixel = 0;
    float result_energy = 0.0f;

    for (int i = 0; i < pixel_sizes.size(); i++) {

        cout << "HINT: pixel scale factor:" << pixel_scale_factors[i] << endl;

        // get max_rms and max_rms_per_area
        vector<vector<float>> result;
        fstream input_file(input_file_path);
        if (!input_file.good()) {
            cerr << "File " << input_file_path << " opened fail." << endl;
            return;
        }
        result = getDataMatrix(input_file);

        // 得到放缩后的matrix
        result = getSubSampledMatrix(result, pixel_scale_factors[i], pixel_scale_factors[i]);

        int r = 0;
        string str;
        max_dif = -FLT_MAX;
        sum_dif = 0.0f;
        rms_dif = 0.0f;
        ave_dif = 0.0f;
        rms_per_area = 0.0f;
        R2 = 0.0f;
        total_energy_err = 0.0f;
        abs_energy_err = 0.0f;
        result_energy = 0.0f;
        int result_row = result.size();
        int result_column = result[0].size();
        pixel = ground_truth_row * ground_truth_column;
        int factor_row = ground_truth_row / result_row;
        int factor_col = ground_truth_column / result_column;


        for (int row = 0; row < ground_truth_row; row++) {
            for (int col = 0; col < ground_truth_column; col++) {

                int r = min(result_row - 1, row / factor_row);
                int c = min(result_column - 1, col / factor_col);
                float result_val = result[r][c];
                float ground_truth_val = ground_truth[row][col];
                float dif = abs(result_val - ground_truth_val);

                rms_dif += (dif * dif);
                max_dif = max(max_dif, dif);
                sum_dif += dif;

                abs_energy_err += (receiver_pixel_area * dif);
                result_energy += receiver_pixel_area * result_val;

            }
        }
        rms_dif = sqrt(rms_dif / pixel);

        rms_per_area = rms_dif / sqrt(heliostat_area) * sqrt(receiver_width * receiver_height);

        ave_dif = sum_dif / pixel;
        R2 = getCoefficientOfDetermination(result, ground_truth);


        total_energy_err = abs(result_energy - ground_truth_energy);

        output_file << pixel_scale_factors[i] << "," << pixel_sizes[i] << ", " << max_dif << ", " << sum_dif << ", " << rms_dif << ", " << ave_dif << ", "
            << pixel
            << "," << R2 << "," << rms_per_area << "," << total_energy_err << ","
            << abs_energy_err << endl;
    }
    output_file.close();
}

void split(string str, vector<float>& v, string spacer) {
    v.clear();
    int pos1, pos2;
    int len = spacer.length();     //记录分隔符的长度
    pos1 = 0;
    pos2 = str.find(spacer);
    while (pos2 != string::npos) {
        v.push_back(atof(str.substr(pos1, pos2 - pos1).c_str()));
        pos1 = pos2 + len;
        pos2 = str.find(spacer, pos1);    // 从str的pos1位置开始搜寻spacer
    }
    if (pos1 != str.length()) //分割最后一个部分
        v.push_back(atof(str.substr(pos1).c_str()));
}

void splitStr(string str, vector<string>& v, string spacer) {

    stringstream ss(str);
    v.clear();
    string s;
    while (getline(ss, s, spacer[0])) {
        if (s.back() == '\r') {
            s.pop_back();
        }
        v.push_back(s);
    }

}

void
combineResults(vector<string>& input_file_lists, vector<string> titles, string target_title,
    string output_file_path) {
    int n = input_file_lists.size();
    cout << "n:" << n << endl;
    vector<stringstream> input_ss(n);
    int index = -1;

    for (int i = 0; i < input_file_lists.size(); i++) {
        fstream file(input_file_lists[i], fstream::in);
        input_ss[i] << file.rdbuf();
        file.close();
    }

    fstream output_file(output_file_path, fstream::out);

    string temp_str = input_ss[0].str();
    int row = count(temp_str.begin(), temp_str.end(), '\n');
    vector<string> strs(input_ss.size());

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < input_ss.size(); j++) {
            getline(input_ss[j], strs[j]);
        }
        if (i == 0) {
            vector<string> title_split;
            splitStr(strs[0], title_split, ",");
            for (int j = 0; j < title_split.size(); j++) {
                if (title_split[j] == target_title) {
                    cout << title_split[j] << endl;
                    index = j;
                    break;
                }
            }
            if (index == -1) {
                return;
            } else {
                cout << " index:" << index << " title:" << title_split[index];
            }
            for (int j = 0; j < titles.size() - 1; j++) {
                output_file << titles[j] << ",";
            }
            output_file << titles.back() << endl;
        } else {
            vector<float> values(input_ss.size(), 0.0f);
            for (int j = 0; j < input_ss.size(); j++) {
                vector<float> v;
                split(strs[j], v, ",");
                values[j] = v[index];
            }
            for (int j = 0; j < values.size() - 1; j++) {
                output_file << values[j] << ",";
            }
            output_file << values.back() << endl;
        }
    }

    output_file.close();
}
vector<vector<float>> getDataMatrix(fstream& file) {
    vector<vector<float>> result;
    stringstream ss;
    ss << file.rdbuf();
    file.close();
    string str;
    while (getline(ss, str)) {
        string val_str;
        vector<float> row_data;
        stringstream row_ss(str);
        while (getline(row_ss, val_str, ',')) {
            float val = atof(val_str.c_str());
            row_data.push_back(val);
        }
        result.push_back(row_data);
    }
    return result;
};

vector<vector<float>> getSubSampledMatrix(vector<vector<float>>& origin, int dh, int dw) {

    vector<vector<float>> result;
    int h = origin.size();
    int w = origin[0].size();
    for (int i = 0; i <= h - dh; i += dh) {
        vector<float> row_data;
        for (int j = 0; j <= w - dw; j += dw) {
            float val = 0.0f;
            for (int di = 0; di < dh; di++) {
                for (int dj = 0; dj < dw; dj++) {
                    val += origin[i + di][j + dj];
                }
            }
            val /= (dh * dw);
            row_data.push_back(val);
        }
        result.push_back(row_data);
    }
    return result;
};