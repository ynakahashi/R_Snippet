match.call(get, call("get", "abc", i = FALSE, p = 3))
## -> get(x = "abc", pos = 3, inherits = FALSE)
fun <- function(x, lower = 0, upper = 1) {
   structure((x - lower) / (upper - lower), CALL = match.call())
}
fun(4 * atan(1), u = pi)
match.call(fun, call("fun", 4, lower = 0, upper = 1))


match.call(fun, call("fun", 4, lower = 0, upper = 1))
call("fun", x = 4, lower = 0, upper = 1)
identical(match.call(fun, call("fun", 4, lower = 0, upper = 1)),
          call("fun", x = 4, lower = 0, upper = 1))

get("fun", x = 4, lower = 0, upper = 1)



eval(match.call(fun, call("fun", 4, lower = 0, upper = 1)))
eval(call("fun", x = 4, lower = 0, upper = 1))
