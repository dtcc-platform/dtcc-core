#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <type_traits>

namespace terrain_mesher::core {

template <typename T>
class Raster {
  public:
    Raster() {
        if constexpr (std::is_floating_point_v<T>) {
            no_data_value_ = std::numeric_limits<T>::quiet_NaN();
        } else {
            no_data_value_ = T{};
        }
    }

    Raster(std::size_t width, std::size_t height) : Raster() { allocate(width, height); }

    Raster(Raster&&) noexcept = default;
    Raster& operator=(Raster&&) noexcept = default;

    Raster(const Raster&) = delete;
    Raster& operator=(const Raster&) = delete;

    void allocate(std::size_t width, std::size_t height) {
        width_ = width;
        height_ = height;
        data_ = std::make_unique<T[]>(width * height);
    }

    void set_all(const T& value) {
        std::fill(get_ptr(), get_ptr() + width_ * height_, value);
    }

    std::size_t get_width() const noexcept { return width_; }
    std::size_t get_height() const noexcept { return height_; }

    void set_pos_x(double value) noexcept { pos_x_ = value; }
    void set_pos_y(double value) noexcept { pos_y_ = value; }
    void set_cell_size(double value) noexcept { cell_size_ = value; }
    void set_no_data_value(T value) noexcept { no_data_value_ = value; }

    double get_pos_x() const noexcept { return pos_x_; }
    double get_pos_y() const noexcept { return pos_y_; }
    double get_cell_size() const noexcept { return cell_size_; }
    T get_no_data_value() const noexcept { return no_data_value_; }

    T* get_ptr() const noexcept { return data_.get(); }
    T* get_ptr(std::size_t row) const noexcept { return data_.get() + row * width_; }

    T& value(std::size_t row, std::size_t column) const noexcept {
        return get_ptr(row)[column];
    }

    double col2x(std::size_t column) const noexcept {
        return pos_x_ + (static_cast<double>(column) + 0.5) * cell_size_;
    }

    double row2y(std::size_t row_top_left) const noexcept {
        const auto row_lower_left = height_ - 1 - row_top_left;
        return pos_y_ + (static_cast<double>(row_lower_left) + 0.5) * cell_size_;
    }

    bool empty() const noexcept { return width_ == 0 || height_ == 0; }

  private:
    std::size_t width_ = 0;
    std::size_t height_ = 0;
    double pos_x_ = 0.0;
    double pos_y_ = 0.0;
    double cell_size_ = 1.0;
    T no_data_value_{};
    std::unique_ptr<T[]> data_;
};

using RasterDouble = Raster<double>;

bool is_no_data(double value, double no_data_value) noexcept;
double sample_nearest_valid_avg(const RasterDouble& src,
                                std::size_t row,
                                std::size_t column,
                                int min_averaging_samples = 1);

namespace detail {

inline double safe_get_pixel(const RasterDouble& src,
                             std::int64_t width,
                             std::int64_t height,
                             std::int64_t row,
                             std::int64_t column) {
    if (row < 0 || column < 0 || row >= height || column >= width) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return src.value(static_cast<std::size_t>(row), static_cast<std::size_t>(column));
}

template <std::size_t Size>
inline double average_nan_arr(const std::array<double, Size>& to_average) {
    double sum = 0.0;
    int count = 0;
    for (double value : to_average) {
        if (!std::isnan(value)) {
            sum += value;
            ++count;
        }
    }
    if (count == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return sum / static_cast<double>(count);
}

inline double subsample_raster_3x3(const RasterDouble& src,
                                   double no_data_value,
                                   std::int64_t width,
                                   std::int64_t height,
                                   std::int64_t row,
                                   std::int64_t column) {
    double center_pixel = safe_get_pixel(src, width, height, row, column);
    std::array<double, 4> cross_pixels = {
        safe_get_pixel(src, width, height, row - 1, column),
        safe_get_pixel(src, width, height, row, column - 1),
        safe_get_pixel(src, width, height, row, column + 1),
        safe_get_pixel(src, width, height, row + 1, column),
    };
    std::array<double, 4> diag_pixels = {
        safe_get_pixel(src, width, height, row - 1, column - 1),
        safe_get_pixel(src, width, height, row - 1, column + 1),
        safe_get_pixel(src, width, height, row + 1, column - 1),
        safe_get_pixel(src, width, height, row + 1, column + 1),
    };

    if (center_pixel == no_data_value) {
        center_pixel = std::numeric_limits<double>::quiet_NaN();
    }
    for (double& value : cross_pixels) {
        if (value == no_data_value) {
            value = std::numeric_limits<double>::quiet_NaN();
        }
    }
    for (double& value : diag_pixels) {
        if (value == no_data_value) {
            value = std::numeric_limits<double>::quiet_NaN();
        }
    }

    const double cross_avg = average_nan_arr(cross_pixels);
    const double diag_avg = average_nan_arr(diag_pixels);
    const std::array<double, 6> weighted = {
        center_pixel,
        center_pixel,
        center_pixel,
        cross_avg,
        cross_avg,
        diag_avg,
    };
    return average_nan_arr(weighted);
}

constexpr int kMaxAveragingSamples = 64;

inline double average_samples(const std::array<double, kMaxAveragingSamples>& values, int count) {
    if (count == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double sum = 0.0;
    for (double value : values) {
        sum += value;
    }
    return sum / static_cast<double>(count);
}

}  // namespace detail

inline bool is_no_data(double value, double no_data_value) noexcept {
    return std::isnan(value) || value == no_data_value;
}

inline double sample_nearest_valid_avg(const RasterDouble& src,
                                       std::size_t row,
                                       std::size_t column,
                                       int min_averaging_samples) {
    min_averaging_samples = std::min(min_averaging_samples, detail::kMaxAveragingSamples);

    const std::int64_t row_i = static_cast<std::int64_t>(row);
    const std::int64_t column_i = static_cast<std::int64_t>(column);
    const std::int64_t width = static_cast<std::int64_t>(src.get_width());
    const std::int64_t height = static_cast<std::int64_t>(src.get_height());
    const std::int64_t max_radius =
        static_cast<std::int64_t>(std::sqrt(static_cast<double>(width * width + height * height)));
    const double no_data_value = src.get_no_data_value();

    double z = 0.0;
    if (row_i < height && column_i < width) {
        z = src.value(row, column);
    }
    if (!is_no_data(z, no_data_value)) {
        return z;
    }

    std::array<double, detail::kMaxAveragingSamples> to_average{};
    int average_count = 0;

    auto put_pixel = [&](int x, int y) {
        const std::int64_t dest_row = row_i + y;
        const std::int64_t dest_column = column_i + x;
        double value =
            detail::subsample_raster_3x3(src, no_data_value, width, height, dest_row, dest_column);
        if (!is_no_data(value, no_data_value) && average_count < detail::kMaxAveragingSamples) {
            to_average[average_count] = value;
            ++average_count;
        }
    };

    for (std::int64_t radius = 2;
         radius <= max_radius && average_count < min_averaging_samples;
         ++radius) {
        std::int64_t x = radius - 1;
        std::int64_t y = 0;
        std::int64_t dx = 1;
        std::int64_t dy = 1;
        std::int64_t err = dx - (radius / 2);

        while (x >= y) {
            put_pixel(static_cast<int>(x), static_cast<int>(y));
            put_pixel(static_cast<int>(y), static_cast<int>(x));
            put_pixel(static_cast<int>(-y), static_cast<int>(x));
            put_pixel(static_cast<int>(-x), static_cast<int>(y));
            put_pixel(static_cast<int>(-x), static_cast<int>(-y));
            put_pixel(static_cast<int>(-y), static_cast<int>(-x));
            put_pixel(static_cast<int>(y), static_cast<int>(-x));
            put_pixel(static_cast<int>(x), static_cast<int>(-y));

            if (err <= 0) {
                ++y;
                err += dy;
                dy += 2;
            } else {
                --x;
                dx += 2;
                err += dx - (radius / 2);
            }
        }
    }

    if (average_count == 1) {
        return to_average[0];
    }
    return detail::average_samples(to_average, average_count);
}

}  // namespace terrain_mesher::core
