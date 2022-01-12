#pragma once

#include <fmt/core.h>
#include <fmt/color.h>

namespace Logger{
    namespace LogLevel{
        struct INFO{
            static std::string get(){
                return fmt::format(fmt::emphasis::bold | fg(fmt::color::green), "[INFO]");
            }
        };
        struct WARNING{
            static std::string get(){
                return fmt::format(fmt::emphasis::bold | fg(fmt::color::yellow), "[WARN]");
            }
        };
        struct ERROR{
            static std::string get(){
                return fmt::format(fmt::emphasis::bold | fg(fmt::color::orange), "[ERR] ");
            }
        };
        struct FATAL{
            static std::string get(){
                return fmt::format(fmt::emphasis::bold | fg(fmt::color::red), "[FAT] ");
            }
        };
        struct DEBUG{
            static std::string get(){
                return fmt::format(fmt::emphasis::bold | fg(fmt::color::orange), "[DEB] ");
            }
        };
    };

    namespace{
        namespace detail{
            struct LoggerImpl{
                static const std::string PREFIX;

                template <typename T>
                void Log(T first){
                    fmt::print(fmt::format("{}\n", first));
                }

                template <typename T, typename ...args>
                void Log(T first, args... next){
                    fmt::print(fmt::format("{}", first));
                    Log(next...);
                }
            };
            const std::string LoggerImpl::PREFIX = "[SIM]";
        }

        static detail::LoggerImpl logImpl;
    }

    template <typename LogLevel = LogLevel::INFO, typename ...args>
    void Log(args... next){
        fmt::print( fmt::format( "{}{}: ", detail::LoggerImpl::PREFIX, LogLevel::get() ) );
        logImpl.Log(next...);
    }
}