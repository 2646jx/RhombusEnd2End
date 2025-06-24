
#ifndef RHOMBUS_STATUSOR_H
#define RHOMBUS_STATUSOR_H


#include "status.h"

namespace antchain::util
    {

        // A statusor holds a Status (in the case of an error), or a value T.
        template <typename T>
        class StatusOr
        {
        public:
            inline StatusOr();

            // Builds from a non-OK status. Crashes if an OK status is specified.
            inline StatusOr(const antchain::util::Status& status);  // NOLINT

            // Builds from the specified value.
            inline StatusOr(const T& value);  // NOLINT
            inline StatusOr(T&& value);       // NOLINT

            // Copy constructor.
            inline StatusOr(const StatusOr& other);

            // Move constructor.
            inline StatusOr(StatusOr&& other);

            // Conversion copy constructor, T must be copy constructible from U.
            template <typename U>
            inline StatusOr(const StatusOr<U>& other);

            // Assignment operator.
            inline const StatusOr& operator=(const StatusOr& other);

            // Conversion assignment operator, T must be assignable from U
            template <typename U>
            inline const StatusOr& operator=(const StatusOr<U>& other);

            // Accessors.
            inline const antchain::util::Status& status() const {
                return status_;
            }

            // Shorthand for status().ok().
            inline bool ok() const {
                return status_.ok();
            }


            // Returns value or crashes if ok() is false.
            inline const T& ValueOrDie() const& {
                if (!ok()) {
                    std::cerr << "Attempting to fetch value of non-OK StatusOr\n";
                    std::cerr << status() << std::endl;
                    std::_Exit(1);
                }
                return value_;
            }
            inline T& ValueOrDie() & {
                if (!ok()) {
                    std::cerr << "Attempting to fetch value of non-OK StatusOr\n";
                    std::cerr << status() << std::endl;
                    std::_Exit(1);
                }
                return value_;
            }
            inline const T&& ValueOrDie() const&& {
                if (!ok()) {
                    std::cerr << "Attempting to fetch value of non-OK StatusOr\n";
                    std::cerr << status() << std::endl;
                    std::_Exit(1);
                }
                return std::move(value_);
            }
            inline T&& ValueOrDie() && {
                if (!ok()) {
                    std::cerr << "Attempting to fetch value of non-OK StatusOr\n";
                    std::cerr << status() << std::endl;
                    std::_Exit(1);
                }
                return std::move(value_);
            }

            template <typename U>
            friend class StatusOr;

        private:
            Status status_;
            T value_;
        };

        template <typename T>
        inline StatusOr<T>::StatusOr()
                : status_(antchain::util::error::UNKNOWN, "") {
        }

        template <typename T>
        inline StatusOr<T>::StatusOr(
                const antchain::util::Status& status) : status_(status) {
            if (status.ok()) {
                std::cerr << "::crypto::tink::util::OkStatus() "
                          << "is not a valid argument to StatusOr\n";
                std::_Exit(1);
            }
        }

        template <typename T>
        inline StatusOr<T>::StatusOr(const T& value) : value_(value) {
        }

        template <typename T>
        inline StatusOr<T>::StatusOr(T&& value) : value_(std::move(value)) {
        }

        template <typename T>
        inline StatusOr<T>::StatusOr(const StatusOr& other)
                : status_(other.status_), value_(other.value_) {
        }

        template <typename T>
        inline StatusOr<T>::StatusOr(StatusOr&& other)
                : status_(other.status_), value_(std::move(other.value_)) {
        }

        template <typename T>
        template <typename U>
        inline StatusOr<T>::StatusOr(const StatusOr<U>& other)
                : status_(other.status_), value_(other.value_) {
        }

        template <typename T>
        inline const StatusOr<T>& StatusOr<T>::operator=(const StatusOr& other) {
            status_ = other.status_;
            if (status_.ok()) {
                value_ = other.value_;
            }
            return *this;
        }

        template <typename T>
        template <typename U>
        inline const StatusOr<T>& StatusOr<T>::operator=(const StatusOr<U>& other) {
            status_ = other.status_;
            if (status_.ok()) {
                value_ = other.value_;
            }
            return *this;
        }
    }


#endif //RHOMBUS_STATUSOR_H
