#pragma once

namespace mymath {
	namespace {
		template<class T>
		constexpr const T& min(const T& a, const T& b) {
			return a < b ? a : b;
		}

		template<class T>
		constexpr const T& max(const T& a, const T& b) noexcept {
			return a > b ? a : b;
		}
		 
	}
}