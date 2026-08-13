// core/Error.h - lightweight Result<T> / Error without relying on C++23 std::expected.
#pragma once

#include <optional>
#include <string>
#include <utility>
#include <variant>

namespace biocad {

// A structured error: a category code plus a human-readable message.
struct Error {
    enum class Code {
        Unknown,
        NotFound,
        InvalidArgument,
        Io,
        Parse,
        Unsupported,
        Conflict,
        Internal,
    };

    Code code{Code::Unknown};
    std::string message;

    Error() = default;
    Error(Code c, std::string m) : code(c), message(std::move(m)) {}

    static Error notFound(std::string m)        { return {Code::NotFound, std::move(m)}; }
    static Error invalidArgument(std::string m) { return {Code::InvalidArgument, std::move(m)}; }
    static Error io(std::string m)              { return {Code::Io, std::move(m)}; }
    static Error parse(std::string m)           { return {Code::Parse, std::move(m)}; }
    static Error unsupported(std::string m)     { return {Code::Unsupported, std::move(m)}; }
    static Error internal(std::string m)        { return {Code::Internal, std::move(m)}; }

    [[nodiscard]] const char* codeName() const {
        switch (code) {
            case Code::NotFound:        return "NotFound";
            case Code::InvalidArgument: return "InvalidArgument";
            case Code::Io:              return "Io";
            case Code::Parse:           return "Parse";
            case Code::Unsupported:     return "Unsupported";
            case Code::Conflict:        return "Conflict";
            case Code::Internal:        return "Internal";
            case Code::Unknown:         return "Unknown";
        }
        return "Unknown";
    }
};

// Result<T>: holds either a value of type T or an Error.
template <typename T>
class Result {
public:
    Result(T value) : data_(std::move(value)) {}            // NOLINT(google-explicit-constructor)
    Result(Error err) : data_(std::move(err)) {}            // NOLINT(google-explicit-constructor)

    [[nodiscard]] bool ok() const { return std::holds_alternative<T>(data_); }
    explicit operator bool() const { return ok(); }

    T& value() { return std::get<T>(data_); }
    const T& value() const { return std::get<T>(data_); }

    [[nodiscard]] const Error& error() const { return std::get<Error>(data_); }

    T valueOr(T fallback) const { return ok() ? std::get<T>(data_) : std::move(fallback); }

private:
    std::variant<T, Error> data_;
};

// Result<void>: success or Error, no payload.
template <>
class Result<void> {
public:
    Result() = default;
    Result(Error err) : err_(std::move(err)) {}             // NOLINT(google-explicit-constructor)

    [[nodiscard]] bool ok() const { return !err_.has_value(); }
    explicit operator bool() const { return ok(); }
    [[nodiscard]] const Error& error() const { return *err_; }

    static Result success() { return Result{}; }

private:
    std::optional<Error> err_{};
};

using Status = Result<void>;

}  // namespace biocad
