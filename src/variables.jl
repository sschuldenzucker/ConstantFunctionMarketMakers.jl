"Module exports some variables as _x etc. to avoid conflicts with 'normal' variables."
module Variables

using Symbolics

export _x, _y, _p, _l, _a

const _x = Symbolics.variable("x")
const _y = Symbolics.variable("y")
const _p = Symbolics.variable("p")
const _l = Symbolics.variable("l")
const _a = Symbolics.variable("a")

end
