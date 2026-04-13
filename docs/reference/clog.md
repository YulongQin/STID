# Logging functions for pipeline execution

A comprehensive set of logging functions for tracking pipeline execution
with different log levels and message types. Supports timestamped
messages with various prefixes for different purposes.

## Usage

``` r
clog(message = "", char = ">>>", level = "INFO")

clog_start()

clog_end()

clog_check()

clog_SS(message)

clog_MS(message)

clog_step(message)

clog_loop(message)

clog_normal(message)

clog_title(message)

clog_warn(message)

clog_error(message, call. = TRUE)
```

## Arguments

- message:

  Character, the log message to display

- char:

  Character, prefix character(s) for the log message (default: "\>\>\>")

- level:

  Character, log level - one of "DEBUG", "INFO", "WARNING", "ERROR",
  "CRITICAL"

- call.:

  Logical, whether to include call information in error messages
  (default: TRUE)

## Value

NULL (invisible), messages are printed to console

## Examples

``` r
if (FALSE) { # \dontrun{
clog("This is an info message")
clog_start()
clog_step("Starting analysis step")
clog_loop("Processing sample 1 of 10")
clog_normal("Normal operation message")
clog_warn("This is a warning")
clog_end()
} # }
```
