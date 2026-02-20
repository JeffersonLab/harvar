#include <cctype>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * Convert a string with hexadecimal digits (optionally separated by spaces)
 * into a vector of bytes (stored as char).
 *
 * Examples of accepted input:
 *   "0A1B2C3D"
 *   "0a 1b 2c 3d"
 */
std::vector<char> hexToBytes(const std::string& hexStr)
{
    std::vector<char> bytes;
    bytes.reserve(hexStr.size() / 2);              // rough upper bound

    auto isHex = [](char c) {
        return std::isxdigit(static_cast<unsigned char>(c));
    };

    std::string nibble;   // collects two hex digits
    nibble.reserve(2);

    for (char c : hexStr) {
        if (std::isspace(static_cast<unsigned char>(c)))   // skip spaces/tabs/newlines
            continue;

        if (!isHex(c))
            throw std::invalid_argument("Non‑hex character in input");

        nibble.push_back(c);
        if (nibble.size() == 2) {
            // parse the two‑digit hex number
            char byte = static_cast<char>(std::strtol(nibble.c_str(), nullptr, 16));
            bytes.push_back(byte);
            nibble.clear();
        }
    }

    if (!nibble.empty())
        throw std::invalid_argument("Odd number of hex digits");

    return bytes;
}

int main()
{
    const std::string hex = "0A 1B 2C 3D";
    try {
        std::vector<char> data = hexToBytes(hex);

        std::cout << "Converted " << data.size() << " bytes:\n";
        for (unsigned char b : data)
            std::cout << std::hex << std::uppercase << static_cast<int>(b) << ' ';
        std::cout << '\n';

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return EXIT_FAILURE;
    }
}
