/*
 * Copyright 2021 Xilinx, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef _XF_SOLVER_MATRIX_UTILS_HPP_
#define _XF_SOLVER_MATRIX_UTILS_HPP_

#include "ap_fixed.h"
#include "utils/x_hls_utils.h"
#include "hls_x_complex.h"
#include "utils/std_complex_utils.h"
#include "hls_math.h"

#ifndef __SYNTHESIS__
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <stdio.h>
#endif

namespace xf {
namespace solver {

struct NoTranspose {
    const static int TransposeType = 0;
    template <int RowsA, int ColsA, typename InputType>
    static InputType GetElement(const InputType A[RowsA][ColsA], unsigned Row, unsigned Col) {
        return A[Row][Col];
    }
};

struct Transpose {
    const static int TransposeType = 1;
    template <int RowsA, int ColsA, typename InputType>
    static InputType GetElement(const InputType A[RowsA][ColsA], unsigned Row, unsigned Col) {
        return A[Col][Row];
    }
};

struct ConjugateTranspose {
    const static int TransposeType = 2;
    template <int RowsA, int ColsA, typename InputType>
    static InputType GetElement(const InputType A[RowsA][ColsA], unsigned Row, unsigned Col) {
        return hls::x_conj(A[Col][Row]);
    }
};

template <class TransposeForm, int RowsA, int ColsA, typename InputType>
InputType GetMatrixElement(const InputType A[RowsA][ColsA], unsigned Row, unsigned Col) {
    return TransposeForm::template GetElement<RowsA, ColsA, InputType>(A, Row, Col);
}

static half x_sqrt(half x) {
    return hls::half_sqrt(x);
}
static float x_sqrt(float x) {
    return sqrtf(x);
}
static double x_sqrt(double x) {
    return sqrt(x);
}
template <int W1, int I1, ap_q_mode Q1, ap_o_mode O1, int N1>
ap_fixed<W1, I1, Q1, O1, N1> x_sqrt(ap_fixed<W1, I1, Q1, O1, N1> x) {
#pragma HLS INLINE
    return hls::sqrt((float)x);
}

static float x_copysign(float a, float b) {
    return copysignf(a, b);
}
static double x_copysign(double a, double b) {
    return copysign(a, b);
}

static float x_sign(float x) {
    if (x == 0.0f) {
        return 0.0f;
    } else {
        return copysignf(1.0f, x);
    }
}
static double x_sign(double x) {
    if (x == 0.0) {
        return 0.0;
    } else {
        return copysign(1.0, x);
    }
}
template <int W1, int I1, ap_q_mode Q1, ap_o_mode O1, int N1>
ap_fixed<W1, I1, Q1, O1, N1> x_sign(ap_fixed<W1, I1, Q1, O1, N1> x) {
    ap_fixed<W1, I1, Q1, O1, N1> tmp = 0;
    if (x > 0) {
        tmp = 1;
    } else {
        tmp = -1;
    }
    return tmp;
}

static half x_rsqrt(half x) {
    return hls::half_rsqrt(x);
}
static float x_rsqrt(float x) {
    return hls::rsqrtf(x);
}
static double x_rsqrt(double x) {
    return hls::rsqrt(x);
}
template <int W1, int I1, ap_q_mode Q1, ap_o_mode O1, int N1>
ap_fixed<W1, I1, Q1, O1, N1> x_rsqrt(ap_fixed<W1, I1, Q1, O1, N1> x) {
#pragma HLS INLINE
    return hls::rsqrt((float)x);
}

static int x_isneg(float x) {
    fp_struct<float> fs = x;
    return fs.__signbit();
}
static int x_isneg(double x) {
    fp_struct<double> fs = x;
    return fs.__signbit();
}

#ifndef __SYNTHESIS__
template <typename T>
struct xil_printer {
    static std::string to_s(T x, unsigned prec = 10) {
        std::stringstream ss;
        ss << std::setiosflags(std::ios::fixed) << std::setprecision(prec) << x;

        return ss.str();
    }
};

template <typename T>
struct xil_printer<hls::x_complex<T> > {
    static std::string to_s(hls::x_complex<T> x, unsigned prec = 10) {
        typedef xil_printer<T> printer;

        std::stringstream ss;
        bool neg_imag = x.imag() <= -0 ? true : false;

        T imag = neg_imag ? (T)-x.imag() : (T)x.imag();

        ss << printer::to_s(x.real(), prec) << (neg_imag ? " - j*" : " + j*") << printer::to_s(imag, prec);
        return ss.str();
    }
};

template <unsigned ROWS, unsigned COLS, typename T, class TransposeForm>
void print_matrix(T a[ROWS][COLS], std::string prefix = "", unsigned prec = 10, unsigned matlab_format = 0) {
    typedef xil_printer<T> printer;

    std::string res[ROWS][COLS];
    unsigned widest_entry = 0;

    for (int r = 0; r < ROWS; r++) {
        for (int c = 0; c < COLS; c++) {
            T tmp;

            tmp = GetMatrixElement<TransposeForm, ROWS, COLS, T>(a, r, c);

            res[r][c] = printer::to_s(tmp, prec);
            if (res[r][c].length() > widest_entry) {
                widest_entry = res[r][c].length();
            }
        }
    }

    char col_gap_str[5] = "    ";
    unsigned col_width = widest_entry;

    for (int r = 0; r < ROWS; r++) {
        if (!matlab_format) {
            printf("%s|", prefix.c_str());
        } else {
            if (r == 0) {
                printf("%s[", prefix.c_str());
            } else {
                printf("%s ", prefix.c_str());
            }
        }
        for (int c = 0; c < COLS; c++) {
            unsigned num_spaces_needed = col_width - res[r][c].length();
            for (int x = 0; x < num_spaces_needed; x++) {
                printf(" ");
            }
            if (!matlab_format) {
                printf("(%s)", res[r][c].c_str());
            } else {
                printf("%s", res[r][c].c_str());
            }

            if (c != COLS - 1) {
                printf("%s", col_gap_str);
                if (matlab_format) {
                    printf(",");
                }
            }
        }
        if (!matlab_format) {
            printf(" |\n");
        } else {
            if (r == ROWS - 1) {
                printf(" ];\n");
            } else {
                printf(" ;\n");
            }
        }
    }
}
#endif

} // namespace solver
} // namespace xf

#endif

