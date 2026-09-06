//
// Created by cd on 6/9/26.
//

#ifndef CORSET_BINARYIO_H
#define CORSET_BINARYIO_H
// src/BinaryIO.h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>
#include <stdexcept>

namespace BinaryIO {

// ── Magic number and version ─────────────────────────────────────────────────
// Written at the start of every file so we can detect corruption or
// version mismatches on load.
static constexpr uint32_t MAGIC   = 0x43525354; // "CRST"
static constexpr uint32_t VERSION = 1;

inline void write_header(std::ostream &out) {
    uint32_t m = MAGIC, v = VERSION;
    out.write(reinterpret_cast<const char*>(&m), sizeof(m));
    out.write(reinterpret_cast<const char*>(&v), sizeof(v));
}

inline void read_header(std::istream &in) {
    uint32_t m, v;
    in.read(reinterpret_cast<char*>(&m), sizeof(m));
    in.read(reinterpret_cast<char*>(&v), sizeof(v));
    if (m != MAGIC)
        throw std::runtime_error("BinaryIO: bad magic — file is corrupt or wrong type");
    if (v != VERSION)
        throw std::runtime_error("BinaryIO: version mismatch");
}

// ── POD scalar ───────────────────────────────────────────────────────────────
template<typename T>
void write_pod(std::ostream &out, const T &val) {
    out.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

template<typename T>
void read_pod(std::istream &in, T &val) {
    in.read(reinterpret_cast<char*>(&val), sizeof(T));
}

// ── string ───────────────────────────────────────────────────────────────────
// Format: [uint32_t length][char * length]
inline void write_string(std::ostream &out, const std::string &s) {
    uint32_t len = static_cast<uint32_t>(s.size());
    out.write(reinterpret_cast<const char*>(&len), sizeof(len));
    out.write(s.data(), len);
}

inline void read_string(std::istream &in, std::string &s) {
    uint32_t len;
    in.read(reinterpret_cast<char*>(&len), sizeof(len));
    s.resize(len);
    in.read(s.data(), len);
}

// ── vector<POD> ──────────────────────────────────────────────────────────────
// Format: [uint64_t count][T * count]
// Only valid for types where sizeof(T) bytes is the complete value —
// i.e. no pointers, no shared_ptr, no std::string members.
template<typename T>
void write_pod_vector(std::ostream &out, const std::vector<T> &v) {
    uint64_t n = v.size();
    out.write(reinterpret_cast<const char*>(&n), sizeof(n));
    if (n > 0)
        out.write(reinterpret_cast<const char*>(v.data()), n * sizeof(T));
}

template<typename T>
void read_pod_vector(std::istream &in, std::vector<T> &v) {
    uint64_t n;
    in.read(reinterpret_cast<char*>(&n), sizeof(n));
    v.resize(n);
    if (n > 0)
        in.read(reinterpret_cast<char*>(v.data()), n * sizeof(T));
}

// ── vector<string> ───────────────────────────────────────────────────────────
inline void write_string_vector(std::ostream &out, const std::vector<std::string> &v) {
    uint64_t n = v.size();
    out.write(reinterpret_cast<const char*>(&n), sizeof(n));
    for (const auto &s : v)
        write_string(out, s);
}

inline void read_string_vector(std::istream &in, std::vector<std::string> &v) {
    uint64_t n;
    in.read(reinterpret_cast<char*>(&n), sizeof(n));
    v.resize(n);
    for (auto &s : v)
        read_string(in, s);
}

// ── vector<vector<T>> (2D, POD inner) ────────────────────────────────────────
template<typename T>
void write_pod_vector2d(std::ostream &out, const std::vector<std::vector<T>> &v) {
    uint64_t n = v.size();
    out.write(reinterpret_cast<const char*>(&n), sizeof(n));
    for (const auto &inner : v)
        write_pod_vector(out, inner);
}

template<typename T>
void read_pod_vector2d(std::istream &in, std::vector<std::vector<T>> &v) {
    uint64_t n;
    in.read(reinterpret_cast<char*>(&n), sizeof(n));
    v.resize(n);
    for (auto &inner : v)
        read_pod_vector(in, inner);
}

// ── vector<vector<vector<T>>> (3D, POD inner) ─────────────────────────────────
template<typename T>
void write_pod_vector3d(std::ostream &out, const std::vector<std::vector<std::vector<T>>> &v) {
    uint64_t n = v.size();
    out.write(reinterpret_cast<const char*>(&n), sizeof(n));
    for (const auto &mid : v)
        write_pod_vector2d(out, mid);
}

template<typename T>
void read_pod_vector3d(std::istream &in, std::vector<std::vector<std::vector<T>>> &v) {
    uint64_t n;
    in.read(reinterpret_cast<char*>(&n), sizeof(n));
    v.resize(n);
    for (auto &mid : v)
        read_pod_vector2d(in, mid);
}

} // namespace BinaryIO
#endif //CORSET_BINARYIO_H
