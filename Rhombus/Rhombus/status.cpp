
#include <sstream>
#include "status.h"

namespace antchain::util
    {
        using namespace std;
        namespace {
            const Status& GetInternal() {
                static const Status status(antchain::util::error::INTERNAL, "");
                return status;
            }

            const Status& GetUnknown() {
                static const Status status(antchain::util::error::UNKNOWN, "");
                return status;
            }

            const Status& GetOk() {
                static const Status status;
                return status;
            }

        }  // namespace

        Status::Status() : code_(antchain::util::error::OK), message_("") { }


        Status::Status(antchain::util::error::Code error, const std::string& error_message)
                : code_(error), message_(error_message)
        {
            if (code_ == antchain::util::error::OK) {
                message_.clear();
            }
        }

//        Status& Status::operator=(const Status& other) {
//            code_ = other.code_;
//            message_ = other.message_;
//            return *this;
//        }


        const Status& Status::INTERNAL = GetInternal();
        const Status& Status::UNKNOWN = GetUnknown();
        const Status& Status::OK = GetOk();


        std::string Status::ToString() const {
            if (code_ == antchain::util::error::OK) {
                return "OK";
            }

            std::ostringstream oss;
            oss << code_ << ": " << message_;
            return oss.str();
        }


        std::string ErrorCodeString(antchain::util::error::Code error) {
            const char *str = "";
            switch (error) {
                case antchain::util::error::OK:
                    str = "OK";
                    break;
                case antchain::util::error::NULL_PTR:
                    str = "NULL_POINTER";
                    break;
                case antchain::util::error::UNKNOWN:
                    str = "UNKNOWN";
                    break;
                case antchain::util::error::INVALID_ARGUMENT:
                    str = "INVALID_ARGUMENT";
                    break;
                case antchain::util::error::UNSUPPORTED_BITS_OR_POLY_MODULUS_DEGREE:
                    str = "UNSUPPORTED_BITS_OR_POLY_MODULUS_DEGREE";
                    break;
                case antchain::util::error::INVALID_VECTOR_SIZE:
                    str = "INVALID_VECTOR_SIZE";
                    break;
                case antchain::util::error::INCORRECT_PUBLIC_KEY_NUMS:
                    str = "INCORRECT_PUBLIC_KEY_NUMS";
                    break;
                case antchain::util::error::INCORRECT_GALOIS_KEY_NUMS:
                    str = "INCORRECT_GALOIS_KEY_NUMS";
                    break;
                case antchain::util::error::INCORRECT_CIPHERTEXT_NUMS:
                    str = "INCORRECT_CIPHERTEXT_NUMS";
                    break;
                case antchain::util::error::INVALID_MATRIX_ROWS_COLS:
                    str = "INVALID_MATRIX_ROWS_COLS";
                    break;
                case antchain::util::error::INVALID_RANDOM_VECTOR_BITS:
                    str = "INVALID_RANDOM_VECTOR_BITS";
                    break;
                case antchain::util::error::UNINITIALIZED_KEY:
                    str = "UNINITIALIZED_KEY";
                    break;
                case antchain::util::error::INTERNAL:
                    str = "INTERNAL";
                    break;
                case antchain::util::error::INVALID_CRT_PARAMS:
                    str = "INVALID_CRT_PARAMS";
                    break;
                case antchain::util::error::UNSUPPORTED_POLY_MODULUS_DEGREE:
                    str = "UNSUPPORTED_POLY_MODULUS_DEGREE";
                    break;
            }
            // Avoid using a "default" in the switch, so that the compiler can
            // give us a warning, but still provide a fallback here.
            return std::string(str);
        }

        extern ostream& operator<<(ostream& os, antchain::util::error::Code code) {
            os << ErrorCodeString(code);
            return os;
        }

        extern ostream& operator<<(ostream& os, const Status& other) {
            os << other.ToString();
            return os;
        }
    }
