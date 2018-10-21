# test_utils.jl : utilities for testing SimpleMultigrid

# log resnorm
function log(mg,ϵ,β)
    println(string("+","-"^30,"+"))
    println(string("|"," "^4,@sprintf("ϵ = %6.3e, β = %2.0f ",ϵ,β)," "^4,"|"))
    println(string("+","-"^9,"+","-"^20,"+"))
    println(string("| V-cycle |  |rʰ|"," "^7,"ratio","  |"))
    println(string("+","-"^9,"+","-"^20,"+"))
    println(string(@sprintf("|   %4i  |",0),@sprintf("  %6.3e",mg.resnorm[1])," "^8," |"))
    for i in 2:length(mg.resnorm)
        println(string(@sprintf("|   %4i  |",i-1),@sprintf("  %6.3e",mg.resnorm[i]),@sprintf("  %3.2f",mg.resnorm[i]/mg.resnorm[i-1]),"   |"))
    end
    println(string("+","-"^9,"+","-"^20,"+"))
end

# log resnorm
function log(mg,d)
    println(string("+","-"^30,"+"))
	println(string("|"," "^11,@sprintf("n = %4i",(size(mg.grids[1].A,1))^(1/d)+1)," "^11,"|"))
    println(string("+","-"^9,"+","-"^20,"+"))
    println(string("| V-cycle |  |rʰ|"," "^7,"ratio","  |"))
    println(string("+","-"^9,"+","-"^20,"+"))
    println(string(@sprintf("|   %4i  |",0),@sprintf("  %6.3e",mg.resnorm[1])," "^8," |"))
    for i in 2:length(mg.resnorm)
        println(string(@sprintf("|   %4i  |",i-1),@sprintf("  %6.3e",mg.resnorm[i]),@sprintf("  %3.2f",mg.resnorm[i]/mg.resnorm[i-1]),"   |"))
    end
    println(string("+","-"^9,"+","-"^20,"+"))
end
