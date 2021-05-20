using Plots
# plotlyjs()

logocolors = Colors.JULIA_LOGO_COLORS
leftBound = -3
rightBound = 3
upperBound = 2
lowerBound = -0.8

barWidth = leftBound:0.1:rightBound
barHeight = zeros(size(barWidth))

cartWidth = 1.5
cartHeight = 0.5

startPointPend = 0
pendLength = 1.0
thetaTest = 0
ballRad = 0.15
wheelRad = 0.15
pendWidth = 0.1
sol = solSwing

drawCartPole(sol(0.0)[1],0,pi/6)

numDigits = ndigits(length(sol.t));
fname = "video/cartpoleSwingUp";
i = 1;

for t in sol.t
    drawCartPole(sol(t)[1],0,sol(t)[2]);
    numZeros = numDigits - ndigits(i);
    firstPart = join(["0" for j in 1:numZeros]);
    imagei = firstPart*string(i);
    fn = fname*imagei;
    savefig(fn)
    i = i + 1;
end

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
pendbar(w, h, x, y, θ) = Shape([x,x+w*cos(θ),x+w*cos(θ)-h*sin(θ),x+w*cos(θ)-h*sin(θ)-w*cos(θ)], [y,y+w*sin(θ),y+w*sin(θ)+h*cos(θ),y+w*sin(θ)+h*cos(θ)-w*sin(θ)])
plot(pendbar(0.2, 1.0, 0.0, 0.0, π/4),aspect_ratio=1.0)

function drawCartPole(x,y,θ)
    # plot ground and cart
    plot(barWidth,barHeight,c=:black,linewidth=4.0, xlims=(leftBound-0.1,rightBound+0.1), ylims=(lowerBound,upperBound), axis=nothing, legend=false, colorbar=false, grid=false, showaxis=false,aspect_ratio=1.0)
    plot!(rectangle(cartWidth,cartHeight,x-cartWidth/2,y+2*wheelRad),color=logocolors.red,xlims=(leftBound-0.1,rightBound+0.1), ylims=(lowerBound,upperBound), axis=nothing,legend=false, colorbar=false, grid=false, showaxis=false,aspect_ratio=1.0)


    # pendx = range(0,stop=pendLength,length=10)
    # pendy = range(0,stop=pendLength,length=10)
    # plot!(x .+ pendx.*cos(θ-pi/2),cartHeight/2 .+ 2*wheelRad .+ pendy.*sin(θ-pi/2),color=:black,markerstrokecolor=:black,linewidth=8.0,axis=nothing,legend=false, colorbar=false, grid=false, showaxis=false,aspect_ratio=1.0)

    # plot wheels
    plot!(circleShape(x - cartWidth/2+wheelRad,wheelRad,wheelRad),color=:black,fill=(0, 1.0, logocolors.green), xlims=(leftBound-0.1,rightBound+0.1), ylims=(lowerBound,upperBound), axis=nothing,legend=false, colorbar=false, grid=false, showaxis=false,aspect_ratio=1.0)
    plot!(circleShape(x .+ cartWidth/2-wheelRad,wheelRad,wheelRad),color=:black,fill=(0, 1.0, logocolors.green),axis=nothing,legend=false, colorbar=false, grid=false, showaxis=false,aspect_ratio=1.0)
    # plot!(circleShape(x - cartWidth/2+wheelRad,wheelRad,0.3*wheelRad),color=:black,fill=(0, 1.0, :white),axis=nothing,legend=false, colorbar=false, grid=false, showaxis=false,aspect_ratio=1.0)
    # plot!(circleShape(x .+ cartWidth/2-wheelRad,wheelRad,0.3*wheelRad),color=:black,fill=(0, 1.0, :white),axis=nothing,legend=false, colorbar=false, grid=false, showaxis=false,aspect_ratio=1.0)

    # plot # plot pendulum and ball at end of pendulum
    plot!(pendbar(pendWidth,pendLength,x-(pendWidth/2)*cos(θ+pi),y+cartHeight/2+2*wheelRad-(pendWidth/2)*sin(θ+pi),θ+pi),color=logocolors.blue,xlims=(leftBound-0.1, rightBound+0.1), ylims=(lowerBound,upperBound),axis=nothing,legend=false, colorbar=false, grid=false, showaxis=false,aspect_ratio=1.0)
    plot!(circleShape(x .+ pendLength*cos(θ-pi/2), cartHeight/2 .+ 2*wheelRad .+ pendLength*sin(θ-pi/2),ballRad),color=:black,fill=(0, 1.0, logocolors.purple), xlims=(leftBound-0.1,rightBound+0.1), ylims=(lowerBound,upperBound), axis=nothing,legend=false, colorbar=false, grid=false, showaxis=false,aspect_ratio=1.0)
    plot!(circleShape(x, cartHeight/2 .+ 2*wheelRad,0.5*ballRad),color=:black,fill=(0, 1.0, :white), xlims=(leftBound-0.1,rightBound+0.1), ylims=(lowerBound,upperBound), axis=nothing,legend=false, colorbar=false, grid=false, showaxis=false,aspect_ratio=1.0)

end

function circleShape(x,y,r)
    th = range(0,stop=2*pi,length=100);
    zetac = x + y*im;
    zeta = zetac .+ r*exp.(im.*th);
    return real.(zeta), imag.(zeta)
end
