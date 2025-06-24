

#ifndef RHOMBUS_STATUS_H
#define RHOMBUS_STATUS_H

#include <string>
#include <iostream>

namespace antchain::util
    {
        //Error code
        namespace error {
            enum Code {
                OK = 0,

                NULL_PTR = 1,

                UNKNOWN = 2,

                INVALID_ARGUMENT = 3,

                UNSUPPORTED_BITS_OR_POLY_MODULUS_DEGREE = 4,

                INVALID_VECTOR_SIZE = 5,

                INCORRECT_PUBLIC_KEY_NUMS = 6,

                INCORRECT_GALOIS_KEY_NUMS = 7,

                INCORRECT_CIPHERTEXT_NUMS = 8,

                INVALID_MATRIX_ROWS_COLS = 9,

                INVALID_RANDOM_VECTOR_BITS = 10,

                UNINITIALIZED_KEY = 11,

                INTERNAL = 13,

                INVALID_CRT_PARAMS = 20,

                UNSUPPORTED_POLY_MODULUS_DEGREE = 21,
            };
        } //namespace error

        class Status
        {
        public:
            // Creates an OK status
            Status();

            // Make a Status from the specified error and message.
            Status(antchain::util::error::Code error,
                   const std::string& error_message);

//            Status& operator=(const Status& other);

            // Some pre-defined Status objects
            static const Status& OK;  // Identical to 0-arg constructor
            static const Status& INTERNAL;
            static const Status& UNKNOWN;

            // Accessors
            bool ok() const {
                return code_ == error::OK;
            }
            int error_code() const {
                return code_;
            }
            antchain::util::error::Code CanonicalCode() const {
                return code_;
            }
            const std::string& error_message() const { return message_; }

            bool operator==(const Status& x) const;
            bool operator!=(const Status& x) const;

            // NoOp
            void IgnoreError() const {
            }

            std::string ToString() const;

        private:
            antchain::util::error::Code code_;
            std::string message_;
        };


        inline bool Status::operator==(const Status& other) const {
            return (this->code_ == other.code_) && (this->message_ == other.message_);
        }

        inline bool Status::operator!=(const Status& other) const {
            return !(*this == other);
        }

        extern std::string ErrorCodeString(antchain::util::error::Code error);

        extern ::std::ostream& operator<<(::std::ostream& os,
                                          antchain::util::error::Code code);
        extern ::std::ostream& operator<<(::std::ostream& os, const Status& other);

        // Returns an OK status, equivalent to a default constructed instance.
        inline Status OkStatus() { return Status(); }

        #define CATCH_EXCEPTION()   \
        catch (std::invalid_argument &e) { return Status(util::error::INVALID_ARGUMENT, std::string(e.what()));   } \
        catch (std::logic_error &e)      { return Status(util::error::INTERNAL, std::string(e.what()));    }        \
        catch (std::runtime_error &e)    { return Status(util::error::INTERNAL, std::string(e.what())); }           \
        catch (...)                 { return Status(util::error::UNKNOWN, "unknown exception."); }

    }



#endif //RHOMBUS_STATUS_H
