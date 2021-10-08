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
    // cout << "test" << endl;
    // test.open("F:/Data/Normal_1/Dis_100/before_smooth_d_2048_r_128.csv", ios::in);

    float origin_pixel = 0.01;
    // 像素缩放系数 2,3,4,5,...,20,x
    vector<int> pixel_scale_factor;
    // 像素大小 0.02,0.03,...,0.20,0.01*x
    vector<float> pixel_sizes;
    // 初始化 像素设置
    // 2000*2000 -> 20*20 ->100
    for (int i = 2; i <= 100; i += 2) {
        pixel_scale_factor.push_back(i);
        pixel_sizes.push_back(origin_pixel * i);
    }
    // 场景参数
    float receiver_height = 20.0;
    float receiver_width = 20.0;
    float heliostat_area = 4.0f * 3.2f;

    // 定日镜 距离 参数
    vector<int > distance = { 20 };
    // for (int i = 20; i <= 100; i += 10) {
    //     distance.push_back(i);
    // }
    // distance.push_back(120);
    // distance.push_back(130);
    // distance.push_back(150);
    // distance.push_back(160);
    // distance.push_back(170);
    // distance.push_back(190);
    // for (int i = 200; i <= 1200; i += 100) {
    //     distance.push_back(i);
    // }

    // 镜面法向
    int normal = 1;

    // 输出文件路径
    string output_file_path_pre = "../Output/Normal_" + to_string(normal) + "/out_";
    string output_file_path_suf = ".csv";

    for (int dis : distance) {
        string input_file_path_pre = "/home/sth/CLionProjects/SolarEnergy_Chier/OutputFiles/GroundTruth_6282/Altitude_90/Normal_" + to_string(normal) + "/Dis_" + to_string(dis) + "/";
        // test.open("F:/Data/Normal_1/Dis_100/before_smooth_d_2048_r_128.csv", ios::in);

        input_file_path_pre = "F:/Data/Normal_" + to_string(normal) + "/Dis_" + to_string(dis) + "/";
        string input_file_path_suf = "before_smooth_d_2048_r_128.csv";
        string output_file_path = output_file_path_pre + to_string(dis) + output_file_path_suf;

        done(origin_pixel, pixel_scale_factor, pixel_sizes, input_file_path_pre + input_file_path_suf,
            output_file_path, heliostat_area, receiver_height, receiver_width);
    }

    string combine_output_path_pre = "../Output/Normal_" + to_string(normal) + "/out_all_";
    string combine_output_path_suf = ".csv";

    string title = "dis_";
    // ----- 1. Calculate the compare result -----
    // 输入文件夹列表，每个文件夹输出一个结果(一个excel表格)
    // // ----- 2. Combine the compare results -----
    vector<string> target = { "rms", "peak_dif", "peak_shift_gt", "peak_shift_c", "rms", "r2", "rms_per_area","total_energy_result","total_energy_err", "total_energy_err_relative", "abs_energy_err" };
    vector<string> output_file_paths;
    for (int dis : distance) {
        output_file_paths.push_back(output_file_path_pre + to_string(dis) + output_file_path_suf);
    }

    string title_pre = title;
    string title_suf = "";
    vector<string> titles;
    for (int dis : distance) {
        titles.push_back(title_pre + to_string(dis) + title_suf);
    }
    cout << "Combine Result:" << endl;
    for (string t : target) {
        string combine_output_path = combine_output_path_pre + t + combine_output_path_suf;
        combineResults(output_file_paths, titles, t, combine_output_path);
    }
    cout << endl;
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
    // 比较结果
    // pixel_scale_factor: 像素缩放因子(2,3,4,5,...,)
    // pixel_size:像素大小(0.02,0.03,...,)
    // peak_dif: 峰值之差 d
    // peak_shift_gt: 峰值偏移,相对于ground truth 数据的peak val(sqrt(delta_x^2 + delta_y^2)) d
    // peak_shift_c: 峰值偏移,相对于理论中心点(接收器中心点) d
    // rms:原始的 rms=sqrt(sum((val_gt-val)^2)/(pixel count)) d
    // ave:平均误差 sum(val_gt-val)/pixel_count d
    // pixel:pixel count d
    // r2: 相关系数 d
    // rms_per_area: 相对rms = sqrt(sum((val_gt-val)^2)*receiver_size / helio_size) d
    // total_energy_gt:
    // total_energy_result:
    // total_energy_err: 总能量之差 sum(val_gt)-sum(val) d
    // total_energy_err_relative: 相对总能量之差 (sum(val_gt)-sum(val))/sum(val_gt) d
    // abs_energy_err:总能量差和(轮廓之间的差体积) pixel_szie*(val_gt-val)
    output_file
        << "pixel_scale_factor,pixel_size,peak_dif,peak_shift_gt,peak_shift_c,rms,ave,pixel,r2,rms_per_area,total_energy_gt,total_energy_result,total_energy_err,total_energy_err_relative,abs_energy_err"
        << endl;


    // output_file
    //     << "pixel_scale_factor,pixel_size,max,sum,rms,ave,pixel,r2,rms_per_area,total_energy_err,abs_energy_err"
    //     << endl;
    cout << "HINT::Input file path:" << input_file_path << endl;

    fstream ground_truth_file(input_file_path, ios::in);
    if (!ground_truth_file.good()) {
        cerr << "Input file " << input_file_path << " opened fail." << endl;
        return;
    }
    vector<vector<float>> ground_truth = getDataMatrix(ground_truth_file);
    ground_truth_file.close();

    // get result data size (h,w)
    int ground_truth_row = ground_truth.size();
    int ground_truth_column = ground_truth[0].size();

    float receiver_pixel_height = receiver_height / ground_truth_row;
    float receiver_pixel_width = receiver_width / ground_truth_column;

    float receiver_pixel_area = receiver_pixel_height * receiver_pixel_width;

    // abs_energy
    float max_dif; // 最大差
    float sum_dif; // 所有像素差的和
    float rms_dif; // 普通rms(与像素个数、定日镜大小 有关)
    float ave_dif; // 平均差(sum_dif/像素个数)
    float R2;       // 相关系数 R^2
    float rms_per_area; // 相对rms(与像素个数、定日镜大小 无关)
    float total_energy_err; // 接收器上总能量的差(与光斑的分布无关，只是单纯的总能量)
    float abs_energy_err; // 所有像素差的和 * 像素的大小
    float total_energy_ground_truth = 0.0f; // 最小像素时的总能量和
    float total_energy_result = 0.0f;   // 像素缩放后的总能量和
    float peak_ground_truth = 0.0;  // 最小像素时的 峰值
    float peak_result = 0.0f;   // 像素缩放后的 峰值
    float peak_diff = 0.0f;     // 峰值的差
    pair<int, int> peak_pos_ground_truth = pair<int, int>(-1, -1);
    pair<int, int> peak_pos_result = pair<int, int>(-1, -1);
    pair<int, int> peak_pos_theory = pair<int, int>(ground_truth_column / 2, ground_truth_row / 2);
    float peak_shift_gt = 0.0f; // 峰值位置 的 偏移
    float peak_shift_c = 0.0f;  // 像素缩放后
    int pixel = 0;              // 像素个数
    float total_energy_err_relative = 0.0f; // 相对的总能量误差(总能量差/总能量)

    for (int row = 0; row < ground_truth_row; row++) {
        for (int col = 0; col < ground_truth_column; col++) {
            total_energy_ground_truth += (receiver_pixel_area * ground_truth[row][col]);
            if (peak_ground_truth < ground_truth[row][col]) {
                peak_ground_truth = ground_truth[row][col];
                peak_pos_ground_truth = pair<int, int>(col, row);
            }
        }
    }
    fstream input_file(input_file_path, ios::in);
    vector<vector<float>> origin;
    if (!input_file.good()) {
        cerr << "File " << input_file_path << " opened fail." << endl;
        return;
    }
    origin = getDataMatrix(input_file);
    input_file.close();

    for (int i = 0; i < pixel_sizes.size(); i++) {
        // cout << "HINT: pixel scale factor:" << pixel_scale_factors[i] << endl;
        // get max_rms and max_rms_per_area
        vector<vector<float>> result;

        cout << "(" << origin.size() << "*" << origin[0].size() << ")";
        // 得到放缩后的matrix
        result = getSubSampledMatrix(origin, pixel_scale_factors[i], pixel_scale_factors[i]);
        // continue;
        cout << "-- " << "*" << pixel_scale_factors[i] << " -->(" << result.size() << "*" << result[0].size() << ")" << endl;

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
        total_energy_result = 0.0f;
        total_energy_err = 0.0f;
        total_energy_err_relative = 0.0f;
        peak_diff = 0.0f;
        peak_result = 0.0f;

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

                total_energy_result += receiver_pixel_area * result_val;
                if (peak_result < result_val) {
                    peak_result = result_val;
                    peak_pos_result = pair<int, int>(col, row);
                }
            }
        }

        rms_dif = sqrt(rms_dif / pixel);
        rms_per_area = rms_dif / sqrt(heliostat_area) * sqrt(receiver_width * receiver_height);
        ave_dif = sum_dif / pixel;

        R2 = getCoefficientOfDetermination(result, ground_truth);
        total_energy_err = abs(total_energy_result - total_energy_ground_truth);
        total_energy_err_relative = total_energy_err / total_energy_ground_truth;
        peak_diff = abs(peak_result - peak_ground_truth);
        // cout << "(" << result.size() << "*" << result[0].size() << ")" << endl;
        // cout << peak_ground_truth << "<--->" << peak_result << " ->" << abs(peak_diff) << " rms_p:" << rms_per_area << endl;

        float dx = (peak_pos_result.first - peak_pos_ground_truth.first) * receiver_pixel_width;
        float dy = (peak_pos_result.second - peak_pos_ground_truth.second) * receiver_pixel_height;

        peak_shift_gt = sqrt(dx * dx + dy * dy);
        dx = (peak_pos_result.first - peak_pos_theory.first) * receiver_pixel_width;
        dy = (peak_pos_result.second - peak_pos_theory.second) * receiver_pixel_height;
        peak_shift_c = sqrt(dx * dx + dy * dy);

        output_file << pixel_scale_factors[i] << "," << pixel_sizes[i] << "," << peak_diff << "," << peak_shift_gt << "," << peak_shift_c << "," << rms_dif
            << "," << ave_dif << "," << pixel << "," << R2 << "," << rms_per_area << "," << total_energy_ground_truth << "," << total_energy_result << "," << total_energy_err << "," << total_energy_err_relative << "," << abs_energy_err << endl;
        // output_file << pixel_scale_factors[i] << "," << pixel_sizes[i] << ", " << max_dif << ", " << sum_dif << ", " << rms_dif << ", " << ave_dif << ", "
        //     << pixel
        //     << "," << R2 << "," << rms_per_area << "," << total_energy_err << ","
        //     << abs_energy_err << endl;
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

                    index = j;
                    break;
                }
            }
            if (index == -1) {
                return;
            } else {
                cout << title_split[index] << ", ";
            }
            output_file << "pixel,";
            for (int j = 0; j < titles.size() - 1; j++) {
                output_file << titles[j] << ",";
            }
            output_file << titles.back() << endl;
        } else {
            vector<float> values(input_ss.size(), 0.0f);
            vector<float> v;
            for (int j = 0; j < input_ss.size(); j++) {

                split(strs[j], v, ",");
                values[j] = v[index];
            }
            output_file << to_string(v[0] * 0.01f).substr(0, 4) << ",";
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
    // file.close();
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