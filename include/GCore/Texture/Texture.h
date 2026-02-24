#pragma once

#include <GCore/api.h>

#include <cmath>
#include <functional>
#include <glm/glm.hpp>
#include <memory>
#include <string>
#include <vector>

RUZINO_NAMESPACE_OPEN_SCOPE

enum class TextureDimension { Unknown = 0, Tex2D = 2, Tex3D = 3 };

class Texture {
   public:
    virtual ~Texture() = default;

    virtual std::shared_ptr<Texture> clone() const = 0;

    virtual TextureDimension get_dimension() const = 0;

    virtual float sample_scalar_2d(const glm::vec2& coord) const
    {
        return 0.0f;
    }
    virtual glm::vec3 sample_vector_2d(const glm::vec2& coord) const
    {
        return glm::vec3(0.0f);
    }
    virtual glm::vec4 sample_rgba_2d(const glm::vec2& coord) const
    {
        return glm::vec4(0.0f);
    }

    virtual float sample_scalar_3d(const glm::vec3& coord) const
    {
        return 0.0f;
    }
    virtual glm::vec3 sample_vector_3d(const glm::vec3& coord) const
    {
        return glm::vec3(0.0f);
    }
    virtual glm::vec4 sample_rgba_3d(const glm::vec3& coord) const
    {
        return glm::vec4(0.0f);
    }

    virtual bool is_vector_texture() const
    {
        return false;
    }
    virtual int get_width() const
    {
        return 0;
    }
    virtual int get_height() const
    {
        return 0;
    }
    virtual int get_depth() const
    {
        return 1;
    }
};

using TextureHandle = std::shared_ptr<Texture>;

class FunctionTexture2D : public Texture {
   public:
    using ScalarFunc = std::function<float(const glm::vec2&)>;
    using VectorFunc = std::function<glm::vec3(const glm::vec2&)>;

    FunctionTexture2D()
        : scalar_func_(nullptr),
          vector_func_(nullptr),
          is_vector_(false)
    {
    }

    explicit FunctionTexture2D(ScalarFunc func)
        : scalar_func_(std::move(func)),
          vector_func_(nullptr),
          is_vector_(false)
    {
    }

    explicit FunctionTexture2D(VectorFunc func)
        : scalar_func_(nullptr),
          vector_func_(std::move(func)),
          is_vector_(true)
    {
    }

    std::shared_ptr<Texture> clone() const override
    {
        auto tex = std::make_shared<FunctionTexture2D>();
        tex->scalar_func_ = scalar_func_;
        tex->vector_func_ = vector_func_;
        tex->is_vector_ = is_vector_;
        return tex;
    }

    TextureDimension get_dimension() const override
    {
        return TextureDimension::Tex2D;
    }

    void set_scalar_function(ScalarFunc func)
    {
        scalar_func_ = std::move(func);
        vector_func_ = nullptr;
        is_vector_ = false;
    }

    void set_vector_function(VectorFunc func)
    {
        scalar_func_ = nullptr;
        vector_func_ = std::move(func);
        is_vector_ = true;
    }

    float sample_scalar_2d(const glm::vec2& coord) const override
    {
        if (scalar_func_)
            return scalar_func_(coord);
        if (vector_func_) {
            glm::vec3 v = vector_func_(coord);
            return (v.x + v.y + v.z) / 3.0f;
        }
        return 0.0f;
    }

    glm::vec3 sample_vector_2d(const glm::vec2& coord) const override
    {
        if (vector_func_)
            return vector_func_(coord);
        if (scalar_func_) {
            float v = scalar_func_(coord);
            return glm::vec3(v, v, v);
        }
        return glm::vec3(0.0f);
    }

    bool is_vector_texture() const override
    {
        return is_vector_;
    }

   private:
    ScalarFunc scalar_func_;
    VectorFunc vector_func_;
    bool is_vector_;
};

class FunctionTexture3D : public Texture {
   public:
    using ScalarFunc = std::function<float(const glm::vec3&)>;
    using VectorFunc = std::function<glm::vec3(const glm::vec3&)>;

    FunctionTexture3D()
        : scalar_func_(nullptr),
          vector_func_(nullptr),
          is_vector_(false)
    {
    }

    explicit FunctionTexture3D(ScalarFunc func)
        : scalar_func_(std::move(func)),
          vector_func_(nullptr),
          is_vector_(false)
    {
    }

    explicit FunctionTexture3D(VectorFunc func)
        : scalar_func_(nullptr),
          vector_func_(std::move(func)),
          is_vector_(true)
    {
    }

    std::shared_ptr<Texture> clone() const override
    {
        auto tex = std::make_shared<FunctionTexture3D>();
        tex->scalar_func_ = scalar_func_;
        tex->vector_func_ = vector_func_;
        tex->is_vector_ = is_vector_;
        return tex;
    }

    TextureDimension get_dimension() const override
    {
        return TextureDimension::Tex3D;
    }

    void set_scalar_function(ScalarFunc func)
    {
        scalar_func_ = std::move(func);
        vector_func_ = nullptr;
        is_vector_ = false;
    }

    void set_vector_function(VectorFunc func)
    {
        scalar_func_ = nullptr;
        vector_func_ = std::move(func);
        is_vector_ = true;
    }

    float sample_scalar_3d(const glm::vec3& coord) const override
    {
        if (scalar_func_)
            return scalar_func_(coord);
        if (vector_func_) {
            glm::vec3 v = vector_func_(coord);
            return (v.x + v.y + v.z) / 3.0f;
        }
        return 0.0f;
    }

    glm::vec3 sample_vector_3d(const glm::vec3& coord) const override
    {
        if (vector_func_)
            return vector_func_(coord);
        if (scalar_func_) {
            float v = scalar_func_(coord);
            return glm::vec3(v, v, v);
        }
        return glm::vec3(0.0f);
    }

    bool is_vector_texture() const override
    {
        return is_vector_;
    }

   private:
    ScalarFunc scalar_func_;
    VectorFunc vector_func_;
    bool is_vector_;
};

class DataTexture2D : public Texture {
   public:
    enum class WrapMode { Clamp, Repeat, Mirror };

    DataTexture2D() : width_(0), height_(0), wrap_mode_(WrapMode::Repeat)
    {
    }

    std::shared_ptr<Texture> clone() const override
    {
        auto tex = std::make_shared<DataTexture2D>();
        tex->width_ = width_;
        tex->height_ = height_;
        tex->wrap_mode_ = wrap_mode_;
        tex->data_ = data_;
        return tex;
    }

    TextureDimension get_dimension() const override
    {
        return TextureDimension::Tex2D;
    }

    void set_wrap_mode(WrapMode mode)
    {
        wrap_mode_ = mode;
    }
    WrapMode get_wrap_mode() const
    {
        return wrap_mode_;
    }

    void resize(int width, int height)
    {
        width_ = width;
        height_ = height;
        data_.resize(width_ * height_ * 4, 0.0f);
    }

    void set_data(const std::vector<float>& data)
    {
        data_ = data;
    }
    const std::vector<float>& get_data() const
    {
        return data_;
    }
    std::vector<float>& get_data_mut()
    {
        return data_;
    }

    void set_pixel(int x, int y, const glm::vec4& value)
    {
        if (x < 0 || x >= width_ || y < 0 || y >= height_)
            return;
        int idx = (y * width_ + x) * 4;
        data_[idx + 0] = value.r;
        data_[idx + 1] = value.g;
        data_[idx + 2] = value.b;
        data_[idx + 3] = value.a;
    }

    glm::vec4 get_pixel(int x, int y) const
    {
        if (x < 0 || x >= width_ || y < 0 || y >= height_)
            return glm::vec4(0.0f);
        int idx = (y * width_ + x) * 4;
        return glm::vec4(
            data_[idx], data_[idx + 1], data_[idx + 2], data_[idx + 3]);
    }

    int get_width() const override
    {
        return width_;
    }
    int get_height() const override
    {
        return height_;
    }

    float sample_scalar_2d(const glm::vec2& coord) const override
    {
        glm::vec4 c = sample_bilinear(coord.x, coord.y);
        return (c.r + c.g + c.b) / 3.0f;
    }

    glm::vec3 sample_vector_2d(const glm::vec2& coord) const override
    {
        glm::vec4 c = sample_bilinear(coord.x, coord.y);
        return glm::vec3(c.r, c.g, c.b);
    }

    glm::vec4 sample_rgba_2d(const glm::vec2& coord) const override
    {
        return sample_bilinear(coord.x, coord.y);
    }

   private:
    glm::vec4 sample_bilinear(float u, float v) const
    {
        u = apply_wrap_mode(u);
        v = apply_wrap_mode(v);

        float fx = u * width_ - 0.5f;
        float fy = v * height_ - 0.5f;

        int x0 = static_cast<int>(std::floor(fx));
        int y0 = static_cast<int>(std::floor(fy));
        int x1 = x0 + 1;
        int y1 = y0 + 1;

        float tx = fx - x0;
        float ty = fy - y0;

        x0 = wrap_coord(x0, width_);
        y0 = wrap_coord(y0, height_);
        x1 = wrap_coord(x1, width_);
        y1 = wrap_coord(y1, height_);

        glm::vec4 c00 = get_pixel(x0, y0);
        glm::vec4 c10 = get_pixel(x1, y0);
        glm::vec4 c01 = get_pixel(x0, y1);
        glm::vec4 c11 = get_pixel(x1, y1);

        glm::vec4 c0 = c00 * (1.0f - tx) + c10 * tx;
        glm::vec4 c1 = c01 * (1.0f - tx) + c11 * tx;

        return c0 * (1.0f - ty) + c1 * ty;
    }

    float apply_wrap_mode(float coord) const
    {
        switch (wrap_mode_) {
            case WrapMode::Clamp: return glm::clamp(coord, 0.0f, 1.0f);
            case WrapMode::Repeat: return coord - std::floor(coord);
            case WrapMode::Mirror: {
                float scaled = coord * 0.5f;
                float integer;
                float frac = std::modf(scaled, &integer);
                int i = static_cast<int>(integer);
                return (i % 2 == 0) ? frac * 2.0f : (1.0f - frac * 2.0f);
            }
        }
        return coord;
    }

    int wrap_coord(int coord, int size) const
    {
        switch (wrap_mode_) {
            case WrapMode::Clamp: return glm::clamp(coord, 0, size - 1);
            case WrapMode::Repeat: return ((coord % size) + size) % size;
            case WrapMode::Mirror: {
                coord = std::abs(coord);
                int period = coord / size;
                return (period % 2 == 0) ? coord % size
                                         : size - 1 - (coord % size);
            }
        }
        return coord;
    }

    std::vector<float> data_;
    int width_;
    int height_;
    WrapMode wrap_mode_;
};

class DataTexture3D : public Texture {
   public:
    enum class WrapMode { Clamp, Repeat, Mirror };

    DataTexture3D()
        : width_(0),
          height_(0),
          depth_(0),
          wrap_mode_(WrapMode::Repeat)
    {
    }

    std::shared_ptr<Texture> clone() const override
    {
        auto tex = std::make_shared<DataTexture3D>();
        tex->width_ = width_;
        tex->height_ = height_;
        tex->depth_ = depth_;
        tex->wrap_mode_ = wrap_mode_;
        tex->data_ = data_;
        return tex;
    }

    TextureDimension get_dimension() const override
    {
        return TextureDimension::Tex3D;
    }

    void set_wrap_mode(WrapMode mode)
    {
        wrap_mode_ = mode;
    }
    WrapMode get_wrap_mode() const
    {
        return wrap_mode_;
    }

    void resize(int width, int height, int depth)
    {
        width_ = width;
        height_ = height;
        depth_ = depth;
        data_.resize(width_ * height_ * depth_ * 4, 0.0f);
    }

    void set_data(const std::vector<float>& data)
    {
        data_ = data;
    }
    const std::vector<float>& get_data() const
    {
        return data_;
    }
    std::vector<float>& get_data_mut()
    {
        return data_;
    }

    void set_pixel(int x, int y, int z, const glm::vec4& value)
    {
        if (x < 0 || x >= width_ || y < 0 || y >= height_ || z < 0 ||
            z >= depth_)
            return;
        int idx = (z * width_ * height_ + y * width_ + x) * 4;
        data_[idx + 0] = value.r;
        data_[idx + 1] = value.g;
        data_[idx + 2] = value.b;
        data_[idx + 3] = value.a;
    }

    glm::vec4 get_pixel(int x, int y, int z) const
    {
        if (x < 0 || x >= width_ || y < 0 || y >= height_ || z < 0 ||
            z >= depth_)
            return glm::vec4(0.0f);
        int idx = (z * width_ * height_ + y * width_ + x) * 4;
        return glm::vec4(
            data_[idx], data_[idx + 1], data_[idx + 2], data_[idx + 3]);
    }

    int get_width() const override
    {
        return width_;
    }
    int get_height() const override
    {
        return height_;
    }
    int get_depth() const override
    {
        return depth_;
    }

    float sample_scalar_3d(const glm::vec3& coord) const override
    {
        glm::vec4 c = sample_trilinear(coord.x, coord.y, coord.z);
        return (c.r + c.g + c.b) / 3.0f;
    }

    glm::vec3 sample_vector_3d(const glm::vec3& coord) const override
    {
        glm::vec4 c = sample_trilinear(coord.x, coord.y, coord.z);
        return glm::vec3(c.r, c.g, c.b);
    }

    glm::vec4 sample_rgba_3d(const glm::vec3& coord) const override
    {
        return sample_trilinear(coord.x, coord.y, coord.z);
    }

   private:
    glm::vec4 sample_trilinear(float u, float v, float w) const
    {
        u = apply_wrap_mode(u);
        v = apply_wrap_mode(v);
        w = apply_wrap_mode(w);

        float fx = u * width_ - 0.5f;
        float fy = v * height_ - 0.5f;
        float fz = w * depth_ - 0.5f;

        int x0 = static_cast<int>(std::floor(fx));
        int y0 = static_cast<int>(std::floor(fy));
        int z0 = static_cast<int>(std::floor(fz));
        int x1 = x0 + 1;
        int y1 = y0 + 1;
        int z1 = z0 + 1;

        float tx = fx - x0;
        float ty = fy - y0;
        float tz = fz - z0;

        x0 = wrap_coord(x0, width_);
        y0 = wrap_coord(y0, height_);
        z0 = wrap_coord(z0, depth_);
        x1 = wrap_coord(x1, width_);
        y1 = wrap_coord(y1, height_);
        z1 = wrap_coord(z1, depth_);

        glm::vec4 c000 = get_pixel(x0, y0, z0);
        glm::vec4 c100 = get_pixel(x1, y0, z0);
        glm::vec4 c010 = get_pixel(x0, y1, z0);
        glm::vec4 c110 = get_pixel(x1, y1, z0);
        glm::vec4 c001 = get_pixel(x0, y0, z1);
        glm::vec4 c101 = get_pixel(x1, y0, z1);
        glm::vec4 c011 = get_pixel(x0, y1, z1);
        glm::vec4 c111 = get_pixel(x1, y1, z1);

        glm::vec4 c00 = c000 * (1.0f - tx) + c100 * tx;
        glm::vec4 c01 = c001 * (1.0f - tx) + c101 * tx;
        glm::vec4 c10 = c010 * (1.0f - tx) + c110 * tx;
        glm::vec4 c11 = c011 * (1.0f - tx) + c111 * tx;

        glm::vec4 c0 = c00 * (1.0f - ty) + c10 * ty;
        glm::vec4 c1 = c01 * (1.0f - ty) + c11 * ty;

        return c0 * (1.0f - tz) + c1 * tz;
    }

    float apply_wrap_mode(float coord) const
    {
        switch (wrap_mode_) {
            case WrapMode::Clamp: return glm::clamp(coord, 0.0f, 1.0f);
            case WrapMode::Repeat: return coord - std::floor(coord);
            case WrapMode::Mirror: {
                float scaled = coord * 0.5f;
                float integer;
                float frac = std::modf(scaled, &integer);
                int i = static_cast<int>(integer);
                return (i % 2 == 0) ? frac * 2.0f : (1.0f - frac * 2.0f);
            }
        }
        return coord;
    }

    int wrap_coord(int coord, int size) const
    {
        switch (wrap_mode_) {
            case WrapMode::Clamp: return glm::clamp(coord, 0, size - 1);
            case WrapMode::Repeat: return ((coord % size) + size) % size;
            case WrapMode::Mirror: {
                coord = std::abs(coord);
                int period = coord / size;
                return (period % 2 == 0) ? coord % size
                                         : size - 1 - (coord % size);
            }
        }
        return coord;
    }

    std::vector<float> data_;
    int width_;
    int height_;
    int depth_;
    WrapMode wrap_mode_;
};

RUZINO_NAMESPACE_CLOSE_SCOPE
