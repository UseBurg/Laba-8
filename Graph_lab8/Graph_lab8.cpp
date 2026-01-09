#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <map>

#define PI 3.1415926535

using namespace std;

struct Pixel {
    unsigned char r, g, b;

    Pixel() : r(0), g(0), b(0) {}
    Pixel(unsigned char r, unsigned char g, unsigned char b) : r(r), g(g), b(b) {}

    Pixel operator*(double scalar) const {
        return Pixel(
            static_cast<unsigned char>(r * scalar),
            static_cast<unsigned char>(g * scalar),
            static_cast<unsigned char>(b * scalar)
        );
    }

    Pixel operator+(const Pixel& other) const {
        return Pixel(
            static_cast<unsigned char>(r + other.r),
            static_cast<unsigned char>(g + other.g),
            static_cast<unsigned char>(b + other.b)
        );
    }

    Pixel& operator+=(const Pixel& other) {
        r = static_cast<unsigned char>(r + other.r);
        g = static_cast<unsigned char>(g + other.g);
        b = static_cast<unsigned char>(b + other.b);
        return *this;
    }
};

struct Image {
    int width;
    int height;
    vector<vector<Pixel>> data;

    Image(int h = 0, int w = 0) : height(h), width(w) {
        if (h > 0 && w > 0) {
            data.resize(h, vector<Pixel>(w));
        }
    }
};

Image read_ppm(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Cannot open file: " + filename);
    }

    string line;
    while (getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] != '#') break;
    }

    if (line != "P3") {
        throw runtime_error("File is not in PPM P3 format.");
    }

    int width, height;
    file >> width >> height;

    int maxVal;
    file >> maxVal;

    Image image(height, width);

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int r, g, b;
            file >> r >> g >> b;
            image.data[y][x] = Pixel(r, g, b);
        }
    }

    cout << "Image '" << filename << "' loaded successfully. Dimensions: "
        << width << "x" << height << "." << endl;

    file.close();
    return image;
}

void write_ppm(const string& filename, const Image& image) {
    ofstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Cannot create file: " + filename);
    }

    file << "P3\n";
    file << image.width << " " << image.height << "\n";
    file << "255\n";

    for (int y = 0; y < image.height; y++) {
        for (int x = 0; x < image.width; x++) {
            const Pixel& p = image.data[y][x];
            file << static_cast<int>(p.r) << " "
                << static_cast<int>(p.g) << " "
                << static_cast<int>(p.b) << " ";
        }
        file << "\n";
    }

    cout << "Image saved to file '" << filename << "'." << endl;
    file.close();
}

Pixel get_pixel_safe(const Image& image, int y, int x) {
    y = max(0, min(y, image.height - 1));
    x = max(0, min(x, image.width - 1));
    return image.data[y][x];
}

Pixel nearest_neighbor(const Image& image, double y, double x) {
    int ny = static_cast<int>(round(y));
    int nx = static_cast<int>(round(x));
    return get_pixel_safe(image, ny, nx);
}

Pixel bilinear(const Image& image, double y, double x) {
    int l = static_cast<int>(floor(x));
    int k = static_cast<int>(floor(y));
    double a = x - l;
    double b = y - k;

    Pixel p1 = get_pixel_safe(image, k, l);
    Pixel p2 = get_pixel_safe(image, k, l + 1);
    Pixel p3 = get_pixel_safe(image, k + 1, l);
    Pixel p4 = get_pixel_safe(image, k + 1, l + 1);

    Pixel result = p1 * ((1 - a) * (1 - b)) +
        p2 * (a * (1 - b)) +
        p3 * ((1 - a) * b) +
        p4 * (a * b);

    return result;
}

double bicubic_kernel(double t) {
    t = fabs(t);
    if (t <= 1) {
        return 1.5 * pow(t, 3) - 2.5 * pow(t, 2) + 1;
    }
    else if (t <= 2) {
        return -0.5 * pow(t, 3) + 2.5 * pow(t, 2) - 4 * t + 2;
    }
    return 0;
}

Pixel bicubic(const Image& image, double y, double x) {
    int l = static_cast<int>(floor(x));
    int k = static_cast<int>(floor(y));

    double a = x - l;
    double b = y - k;

    Pixel result(0, 0, 0);

    for (int m = -1; m <= 2; m++) {
        for (int n = -1; n <= 2; n++) {
            double weight_x = bicubic_kernel(n - a);
            double weight_y = bicubic_kernel(m - b);

            Pixel p = get_pixel_safe(image, k + m, l + n);

            result = result + p * (weight_x * weight_y);
        }
    }

    result.r = max(0, min(255, static_cast<int>(result.r)));
    result.g = max(0, min(255, static_cast<int>(result.g)));
    result.b = max(0, min(255, static_cast<int>(result.b)));

    return result;
}

struct Matrix3x3 {
    double data[3][3];

    Matrix3x3() {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                data[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
    }

    Matrix3x3 operator*(const Matrix3x3& other) const {
        Matrix3x3 result;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                result.data[i][j] = 0;
                for (int k = 0; k < 3; k++) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }
};

Matrix3x3 get_scale_matrix(double sx, double sy) {
    Matrix3x3 m;
    m.data[0][0] = sx;
    m.data[1][1] = sy;
    return m;
}

Matrix3x3 get_rotation_matrix(double angle_deg) {
    double angle_rad = angle_deg * PI / 180.0;
    double cos_a = cos(angle_rad);
    double sin_a = sin(angle_rad);

    Matrix3x3 m;
    m.data[0][0] = cos_a;
    m.data[0][1] = -sin_a;
    m.data[1][0] = sin_a;
    m.data[1][1] = cos_a;
    return m;
}

Matrix3x3 get_shear_matrix(double shx, double shy) {
    Matrix3x3 m;
    m.data[0][1] = shx;
    m.data[1][0] = shy;
    return m;
}

Matrix3x3 inverse_matrix(const Matrix3x3& m) {
    Matrix3x3 inv;

    double det = m.data[0][0] * (m.data[1][1] * m.data[2][2] - m.data[2][1] * m.data[1][2]) -
        m.data[0][1] * (m.data[1][0] * m.data[2][2] - m.data[1][2] * m.data[2][0]) +
        m.data[0][2] * (m.data[1][0] * m.data[2][1] - m.data[1][1] * m.data[2][0]);

    if (fabs(det) < 1e-10) {
        throw runtime_error("Matrix is singular, cannot find inverse.");
    }

    double inv_det = 1.0 / det;

    inv.data[0][0] = (m.data[1][1] * m.data[2][2] - m.data[2][1] * m.data[1][2]) * inv_det;
    inv.data[0][1] = (m.data[0][2] * m.data[2][1] - m.data[0][1] * m.data[2][2]) * inv_det;
    inv.data[0][2] = (m.data[0][1] * m.data[1][2] - m.data[0][2] * m.data[1][1]) * inv_det;

    inv.data[1][0] = (m.data[1][2] * m.data[2][0] - m.data[1][0] * m.data[2][2]) * inv_det;
    inv.data[1][1] = (m.data[0][0] * m.data[2][2] - m.data[0][2] * m.data[2][0]) * inv_det;
    inv.data[1][2] = (m.data[1][0] * m.data[0][2] - m.data[0][0] * m.data[1][2]) * inv_det;

    inv.data[2][0] = (m.data[1][0] * m.data[2][1] - m.data[2][0] * m.data[1][1]) * inv_det;
    inv.data[2][1] = (m.data[2][0] * m.data[0][1] - m.data[0][0] * m.data[2][1]) * inv_det;
    inv.data[2][2] = (m.data[0][0] * m.data[1][1] - m.data[1][0] * m.data[0][1]) * inv_det;

    return inv;
}

void vector_multiply(const Matrix3x3& m, double vec[3], double result[3]) {
    for (int i = 0; i < 3; i++) {
        result[i] = 0;
        for (int j = 0; j < 3; j++) {
            result[i] += vec[j] * m.data[j][i];
        }
    }
}

Image apply_transform(const Image& image, const Matrix3x3& matrix,
    Pixel(*interpolation_fn)(const Image&, double, double)) {

    Matrix3x3 inverse_matrix_val;
    try {
        inverse_matrix_val = inverse_matrix(matrix);
    }
    catch (const runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
        return image;
    }

    double corners[4][3] = {
        {0, 0, 1},
        {image.width - 1, 0, 1},
        {0, image.height - 1, 1},
        {image.width - 1, image.height - 1, 1}
    };

    double min_x = INFINITY, max_x = -INFINITY;
    double min_y = INFINITY, max_y = -INFINITY;

    for (int i = 0; i < 4; i++) {
        double transformed[3];
        vector_multiply(matrix, corners[i], transformed);

        min_x = min(min_x, transformed[0]);
        max_x = max(max_x, transformed[0]);
        min_y = min(min_y, transformed[1]);
        max_y = max(max_y, transformed[1]);
    }

    int new_width = static_cast<int>(ceil(max_x - min_x));
    int new_height = static_cast<int>(ceil(max_y - min_y));

    Image new_image(new_height, new_width);

    double cx = image.width / 2.0;
    double cy = image.height / 2.0;
    double new_cx = new_width / 2.0;
    double new_cy = new_height / 2.0;

    for (int y_new = 0; y_new < new_height; y_new++) {
        for (int x_new = 0; x_new < new_width; x_new++) {
            double vec_new[3] = { x_new - new_cx, y_new - new_cy, 1 };
            double vec_old[3];
            vector_multiply(inverse_matrix_val, vec_new, vec_old);

            double x_old = vec_old[0] + cx;
            double y_old = vec_old[1] + cy;

            if (x_old >= 0 && x_old < image.width && y_old >= 0 && y_old < image.height) {
                new_image.data[y_new][x_new] = interpolation_fn(image, y_old, x_old);
            }
        }
    }

    return new_image;
}

int main() {
    const string INPUT_FILENAME = "image.ppm";
    const string OUTPUT_PREFIX = "result";

    ifstream test_file(INPUT_FILENAME);
    if (!test_file.is_open()) {
        cerr << "Error: Input file '" << INPUT_FILENAME << "' not found." << endl;
        return 1;
    }
    test_file.close();

    Image source_image;
    try {
        source_image = read_ppm(INPUT_FILENAME);
    }
    catch (const exception& e) {
        cerr << "Error reading file: " << e.what() << endl;
        return 1;
    }

    map<string, Pixel(*)(const Image&, double, double)> interpolations = {
        {"near", nearest_neighbor},
        {"bline", bilinear},
        {"bcub", bicubic}
    };

    map<string, Matrix3x3> transformations = {
        {"scale(2.5,0.75)", get_scale_matrix(2.5, 0.75)},
        {"scale-0.5", get_scale_matrix(0.5, 0.5)},
        {"rotation+45", get_rotation_matrix(45)},
        {"rotation-45", get_rotation_matrix(-45)},
        {"shear(0.2,0)", get_shear_matrix(0.2, 0)},
        {"shear(0,0.7)", get_shear_matrix(0, 0.7)},
        {"rotation+45_and_scale-2", get_rotation_matrix(45) * get_scale_matrix(2, 2)}
    };

    for (const auto& t_pair : transformations) {
        for (const auto& i_pair : interpolations) {
            
            Image result_image = apply_transform(source_image, t_pair.second, i_pair.second);

            string output_filename = OUTPUT_PREFIX + "_" + t_pair.first + "_" + i_pair.first + ".ppm";

            try {
                write_ppm(output_filename, result_image);
            }
            catch (const exception& e) {
                cerr << "Error saving file: " << e.what() << endl;
            }
        }
    }

    cout << "All transformations completed." << endl;

    return 0;
}