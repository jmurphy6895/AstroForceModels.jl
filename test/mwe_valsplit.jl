using ValSplit

@valsplit function soundof(Val(animal::Symbol); size::Number=10.0)
    println(size)
    error("Sound not defined for animal: \$animal")
end

function soundof(animal::Val{:dog}; size::Number=10.0)
    if size > 20.0
        return "Woof"
    else
        return "Arf"
    end
end
    
function soundof(animal::Val{:cat}; size=10.0)
    return "Nyan"
end


soundof(:dog)
soundof(:cat)
soundof(:dog; size=30.0)

soundof(:lizard, size=100.0)