///
/// Created by F.Moitzi on 09.01.2022.
///
/// The memory layout of this class is designed to be the same in {\bf Fortran}.
///

#ifndef MULTIARRAY_HPP
#define MULTIARRAY_HPP

#include <cstddef>
#include <type_traits>
#include <algorithm>
#include <numeric>
#include <array>

namespace detail {
  template<bool...>
  struct bool_pack;
  template<bool... bs>
  using all_true = std::is_same<bool_pack<bs..., true>, bool_pack<true, bs...>>;
}

template<typename... Ts>
using all_true = detail::all_true<Ts::value...>;

template<typename T, typename... Ts>
using all_same = all_true<std::is_same<T, Ts>...>;


namespace lsms {


  template<typename std::size_t Size>
  using array_t = std::array<std::size_t, Size>;

  template<typename T>
  using allocator_t = std::allocator<T>;

  using size_t = std::size_t;


  class ColumnMajor {

  public:

    template<size_t N>
    static void calculate_strides(const array_t<N> &shape, array_t<N> &strides) {
      size_t size = 1;
      for (size_t i = 0; i < N; ++i) {
        strides[i] = size;
        size *= shape[i];
      }
    }

  };

  template<typename Pointer>
  inline typename std::pointer_traits<Pointer>::element_type *to_plain_pointer(Pointer ptr) {
    return (std::addressof(*ptr));
  }


  template<size_t N, typename... Indices>
  inline size_t get_offset(array_t<N> strides, Indices... indices) {

    const array_t<N> inds{indices...};

    size_t off{0};
    for (auto i = 0; i < N; i++) {
      off += strides[i] * inds[i];
    }

    return off;
  }


  template<typename T, size_t N>
  class NDArray {

  private:

    T *data_ = nullptr; // null pointer literal

    size_t size_{0}; // total number of elements

    array_t<N> shape_; // shape of the n-dims array (can be change)

    array_t<N> strides_; //

    allocator_t<T> allocator_; // allocator

  public:


    NDArray() {
      std::fill(std::begin(shape_), std::end(shape_), 0);
      std::fill(std::begin(strides_), std::end(strides_), 0);
    };

    template<typename... Args>
    explicit NDArray(size_t i, Args... args) : shape_{i, args...} {
      size_ = std::accumulate(std::begin(shape_), std::end(shape_), 1, std::multiplies<decltype(N)>());
      ColumnMajor::calculate_strides<N>(shape_, strides_);
      data_ = allocator_.allocate(size_);
    };

    ~NDArray() {
      allocator_.deallocate(data_, size_);
    }

    NDArray(const NDArray &other)
        : data_(nullptr),
          size_(other.size()),
          shape_{other.shape_},
          strides_{other.strides_} {
      data_ = allocator_.allocate(size_);
      std::copy(other.data_, other.data_ + other.size_, data_);
    }

    friend void swap(NDArray &first, NDArray &second) {
      std::swap(first.data_, second.data_);
      std::swap(first.size_, second.size_);
      std::swap(first.shape_, second.shape_);
      std::swap(first.strides_, second.strides_);
    }

    NDArray &operator=(NDArray other) {
      swap(*this, other);
      return *this;
    }

    // Performant way of how to implement it
    NDArray(NDArray &&other) noexcept
        : data_(std::exchange(other.data_, nullptr)),
          size_(std::move(other.size())),
          shape_(std::move(other.shape_)),
          strides_(std::move(other.strides_)) {
    };

    NDArray &operator=(NDArray &&other) noexcept {

      if (data_) allocator_.deallocate(data_, size_);

      data_ = std::exchange(other.data_, nullptr);
      size_ = std::move(other.size());
      shape_ = std::move(other.shape_);
      strides_ = std::move(other.strides_);

      return *this;
    };


    // get access operator
    template<typename... Indices>
    inline T &operator()(Indices... indices) {
      return data_[get_offset<N>(strides_, indices...)];
    }

    template<typename... Indices>
    inline const T &operator()(Indices... indices) const {
      return data_[get_offset<N>(strides_, indices...)];
    }

    template<typename... Indices>
    inline T operator()(Indices... indices) const {
      return data_[get_offset<N>(strides_, indices...)];
    }

    template<typename Indices>
    inline T &operator[](Indices i) {
      return data_[i];
    }

    std::size_t size() const {
      return size_;
    }

    const array_t<N> &shape() const { return shape_; }

    const array_t<N> &strides() const { return strides_; }

    T *data() {
      return to_plain_pointer(data_);
    }

    const T *data() const {
      return to_plain_pointer(data_);
    }

    size_t shape(size_t i) const { return shape_[i]; }

  };

}

#endif //MULTIARRAY_HPP
