### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ a2fa5231-07fe-448c-bbcb-a42bf790e7bc
begin
    import Pkg
	  Pkg.activate(@__DIR__)
	using OrdinaryDiffEq, Plots, LaTeXStrings, PlutoUI, Symbolics #alternativamente también se puede usar SymEngine pero hay que definir el jacobiano a mano
end

# ╔═╡ 183a6460-a9ca-11eb-340e-872e08f29a9e
using LinearAlgebra

# ╔═╡ 5c4af6e0-022e-11eb-152f-2153b9b672c1
md"""# Ejemplos de sistemas a controlar 

Primer ejemplo

$$J \ddot \theta + MgLsin(\theta)=u$$

"""

# ╔═╡ 85c03440-0274-11eb-3c81-514259334750
md"## Simulamos y pintamos"

# ╔═╡ e3e9d0a0-022c-11eb-04ba-97dcb9ba908a
md"""# Grado relativo 
La idea es "derivar hasta que sale u" 
"""

# ╔═╡ 6831eb02-023b-11eb-0842-fbe8f8eb8937
md"Vamos a usar cálculo simbólico"

# ╔═╡ 1c0442a0-a993-11eb-1b2d-55eb8b868dc7
#algunas definiciones útiles
begin
	#using Latexify  # si no lo tenéis hay que instalarlo
    #Base.show(io::IO, ::MIME"text/html", x::Symbolics.Num) = print(io, latexify(x))
    #Base.show(io::IO, ::MIME"text/html", x::Symbolics.Symbolic) = print(io, latexify(x))
    #Base.show(io::IO, ::MIME"text/html", x::Symbolics.Vector{Equation}) = print(io, latexify(x))
    #Base.show(io::IO, ::MIME"text/html", x::Symbolics.AbstractArray{Num}) = print(io, latexify(x))
	
	Base.diff(f,variable) = expand_derivatives(Symbolics.derivative(f,variable))
	subs(expr,dict) = substitute(expr, dict)
	expand(x) = simplify(x;expand=true) #En ciertas versiones es expand = true
	#solve_for(ecuacion,variable) = Symbolics.solve_for(ecuacion, variable)
	#Han roto symbolics lo programo a mano
	function solve_for(ecuación,variable)
	    ecuación2=Num(ecuación.lhs-ecuación.rhs)
        A=diff(ecuación2,variable)
	    if(diff(A,variable)!=0)
		    error("La ecuación no es lineal")
	    end
	    B=ecuación2-A*variable
	    solución=simplify(expand(-B/A))
    end
	jacobian(f,x) = Symbolics.jacobian(f,x)
	
	md"Definimos las funiones habituales..."
end



# ╔═╡ 3ff84c30-0248-11eb-33ba-f99304ee8672
begin
	#D=x1=x2=t=u=J=M=g=L=0  #eso es cosa de pluto y las macros (si no no es necesario)
	@variables D x1 x2 t u J M g L  # =symbols("D,x1,x2,t,u,J, M, g, L")
end

# ╔═╡ 6d535e50-a99a-11eb-08c1-4d98301e53b3
begin
	#escribir las ecuaciones simbólicas
	x=[x1;x2]
	dx1=x2
	dx2=-M*g*L/J*sin(x1)+u/J
	dx=[dx1;dx2]
end

# ╔═╡ 5525d1a0-a9c7-11eb-1de0-1700c9c20460
y=x1

# ╔═╡ 65893910-a9c7-11eb-2f9f-ebadf375cdbc
dy=diff(y,x1)*dx1+diff(y,x2)*dx2

# ╔═╡ a0571620-a9c7-11eb-070f-297fa9e406b2
ddy=diff(dy,x1)*dx1+diff(dy,x2)*dx2

# ╔═╡ 39ce8af0-023d-11eb-3178-857433b369ae
md"ya aparece u --> grado relativo 2"

# ╔═╡ 4c478bf0-023d-11eb-0dce-439d408010fa
md"¿Es repetitivo podemos hacerlo en general?
 hagamos una **función**"

# ╔═╡ b7e08f60-023d-11eb-13cb-ad5285e900d8
begin
	function derivada(funcion,variables,derivadas)
		d=0
		for i=1:length(variables)
			d=d+diff(funcion,variables[i])*derivadas[i]
		end
		return d
	end
	
	function derivada(f,x,dx,n) #multiple distpach #esto es nuevo
		d=f
		for i=1:n
			d=derivada(d,x,dx)
		end
		return d
	end
end

# ╔═╡ 697f24d0-a99a-11eb-26e6-a749fb04cb0f
#usar la función
der=derivada(y,x,dx,2)

# ╔═╡ 9ecc5ff0-022d-11eb-3715-51831818d04f
md"""# Cambio de variables
vamos a cambiar u por v para ello despejamos de
$$y^{(\gamma)}=v$$
"""

# ╔═╡ 65ad8360-a99a-11eb-3764-aded1015107e
begin
	#definir v
	#v=0
	@variables v
end

# ╔═╡ 4b753d20-a9c8-11eb-320a-d58ec5e3c431
ecuación=der~v

# ╔═╡ 6e31c180-a9c8-11eb-3893-f55415fa4fad
uc=solve_for(ecuación,u)

# ╔═╡ 478bc500-2c06-4508-8794-0d99fee916df
test=solve_for(ecuación,u)

# ╔═╡ 561e07e2-0250-11eb-1380-4b5dfa861cd6
md"vamos a ver que lo hemos hecho bien"

# ╔═╡ 64db3f40-a99a-11eb-01df-d3ca4846ca33
#sustituir u en la derivada y ver que sale v
expand(subs(der,u=>uc))

# ╔═╡ de497f50-0250-11eb-028e-ffe56cf1ba45
md"La idea es controlar un nuevo sistema:

$$\dot z_1=z_2$$

$$...$$

$$\dot z_{\gamma}=v$$

$$y=z_1$$
"

# ╔═╡ c2f97c40-a99a-11eb-30e7-359effec0a1d
#define z1 y su derivada z2
z1=y

# ╔═╡ 85e1f1f0-a9c9-11eb-2a97-1944f7fcb34b
z2=derivada(y,x,dx)

# ╔═╡ 0db8c9c0-0252-11eb-199c-2d69997c6769
md"veamos que lo hemos hecho bien y la derivada es la corecta (ha de salir v)"

# ╔═╡ dd56c0c2-a99a-11eb-0c87-8f8e2c6321a8
expand(subs(derivada(z2,x,dx), u=>uc))

# ╔═╡ 72e4e740-022d-11eb-3a63-31023b9d464f
md"""# Difeomorfismo
A partir de ahora controlaremos las **z** pero *OJO* que nuestro objetivo son las **x** 

¿Al controlar z **se controla x**?

-Si se puede despejar x(x) entonces si-->difeomorfismo (al menos local)

Para verlo calculamos el jacobiano y vemos si es singular

$$\begin{pmatrix}

\frac{\partial z_1}{\partial x_1}... & \frac{\partial z_1}{\partial x_n} \\
...\\
\frac{\partial z_n}{\partial x_1}... & \frac{\partial z_n}{\partial x_n}

\end{pmatrix}$$

"""

# ╔═╡ f1e0ee80-a99a-11eb-0c06-85f8d7086aac
#definir el vector z, calcular el jacobiano y derivar
z=[z1,z2]

# ╔═╡ 03abf0e0-a9ca-11eb-339c-5599c951871c
jac=jacobian(z,x)

# ╔═╡ 2d9109e0-a9ca-11eb-08a5-7596eee9e4b2
det(jac)

# ╔═╡ 9ec274e0-022d-11eb-085a-7907b8f1eb93
md"""# Estabilización
El sitema en las zs es lineal
Para estabilizar la dinámica se realimentan negativamente las zs con ganancias
$$v=-d_0z_1-d_1z_2-...-d_{\gamma-1}z_{\gamma}$$

Recordemos que $$\dot z_1=z_2$$, $$\ddot z_1= \dot z_2=z_3$$... 
$$z_1^{(\gamma)}=v=-d_0z_1-d_1z_2-...-d_{\gamma-1}z_{\gamma}$$

$$z_1^{(\gamma)} +d_{\gamma-1}z^{\gamma-1}+...+d_1\dot z+ d_0z  =0$$

Que es un sistema lineal con polinomio característico

$$s^{(\gamma)} +d_{\gamma-1}s+...+d_1s+ d_0  =0$$

## Elección de polos

Hay que hacer que el polinomio sea hurwitz
Lo más fácil es constuirlo, si queremos unos polos $p_1$, $p_2$... $p_{\gamma}$ podemos hacer que el polinomio sea 

$$(s-p_1)(s-p_2)...(s-p_{\gamma})$$
"""


# ╔═╡ aed8d350-0265-11eb-1223-1dd0651632b6
begin
	#expandir (s-p1)*(s-p2)) y comparar con lo deseado
	p1=p2=-1
	@variables s
	expand((s-p1)*(s-p2))
end

# ╔═╡ 0d736c40-0266-11eb-3433-ad75c1617314
md"Entonces $d_0$ y $d_1$ son 1 y 2 ;)"

# ╔═╡ 454272d0-a9cb-11eb-36db-ffb885a95a52
string(uc)

# ╔═╡ f54ff5b2-0266-11eb-1df3-49f3c95729c7
md"Ahora lo ponemos todo junto"

# ╔═╡ 2814f040-0267-11eb-3f4d-2933c0c3528a
function control_1_estabilizacion(x,p,t)
	x1,x2=x
	M,g,L,J=p  
	d0=1
	d1=2
	
	y=x1 #o si lo preferís z1
	dy=x2 #z2
	v=-d0*y-d1*dy
	u=J*v + L*M*g*sin(x1)
	return u
end	

# ╔═╡ 3b125190-0269-11eb-2473-33cbb7451779
md"Si queremos podemos sustituir todo **ojo a veces no compensa** (por que queda muy largo y difícil de depurar) en este caso si"

# ╔═╡ 639dce00-0269-11eb-191e-d16a91f3bcf9
u_sustituido=subs(uc,[v=>-2*y-3*dy,L=>1,M=>1,g=>1,J=>1])

# ╔═╡ e0b21ca0-0a2c-11eb-0e50-efa6c355e7c5
md"Pero **ojo** se pierde la interpretación física de lo que hace el control"

# ╔═╡ f162f3f0-022d-11eb-30d9-55fb6e39c8c0
md"""# Seguimiento
Una vez linealizado seguir una referencia es fácil, si definimos errores

$$e_0=y_r-y$$

$$\dot e_0= \dot y_r - \dot y = e_1$$

...

$$e_0^{(\gamma)}=y_r^{(\gamma)} - v$$

Es todo igual pero añadiendo $y_r$ y con el signo contrario (aparece -v en vez de v) es decir:

$$v=d_0e_0+d_1e_1+...+d_{\gamma-1}e_{\gamma-1} + y_r^{(\gamma)}$$

**OJO a la notación** no solo cambian los signos sino que *los índices ahora están "bien"* ya que $y=e_0$ no como antes que era $z_1$
"""

# ╔═╡ a7745f70-026b-11eb-1abe-d1944165d231
function referencia_1_a(t)
	yr=sin(t)
	dyr=cos(t)
	ddyr=-sin(t)
	return [yr,dyr,ddyr]
end

# ╔═╡ a0e6f960-026b-11eb-23c9-410e6d61c4a1
function control_1_seguimiento(x,p,t)
	x1,x2=x
	yr,dyr,ddyr=referencia_1_a(t)
	
	M,g,L,J=p
	
	d0=2
	d1=3
	
	y=x1
	dy=x2
	
	e0=yr-y #esto es lo nuevo
	e1=dyr-dy
	v=d0*e0 + d1*e1 + ddyr #y aquí también aparece la referencia y cambian los signos
	#v=-d0*y-d1*dy esto es lo de antes 
	#v=-d0*z1-d1*z2 o en términos de las zetas
	
	
	u=J*(v + L*M*g*sin(x1)/J)
	return u
end	

# ╔═╡ 506a20d0-026f-11eb-2bde-f94724bd27f3
begin
	#control_1(x,p,t)=0
	#referencia_1(t)=0
	
	#control_1=control_1_estabilizacion
		
	control_1=control_1_seguimiento
	referencia_1=referencia_1_a
end

# ╔═╡ 85c0d350-022e-11eb-0873-c118ffde014c
function derivadas_ejemplo1(x,p,t)
	x1,x2=x
	M,g,L,J=p
	u=control_1(x,p,t)#calculamos el control
	dx1=x2
    dx2=-M*g*L/J*sin(x1)+u/J
	dx=[dx1 dx2] #no hace falta return dx
end

# ╔═╡ 2a98e6e0-a227-11eb-026a-bb5831f344bb
let
	x0 = [0.1 0.1]
	parametros=[1 1 1 1] # una g un poco rara ¿no?
	tspan = (0.0,10.0)
	prob = ODEProblem(derivadas_ejemplo1,x0,tspan,parametros);
	sol = solve(prob,Tsit5());
	
	# Salida y referencia vars=(0,1) significa pintar variable 1 (el x1) vs variable 0 (el tiempo), cambialo para mostrar otras variables si es necesario
	fig1=plot(sol,vars=(0,1), xaxis=L"t",yaxis=L"y(t)",label=L"y(t)")
	#referencia, x1 no se usa para nada pero espera que le pases algún estado
	fig1=plot!(fig1,sol,vars=( (t,x1)->(t,referencia_1(t)[1]) , 0, 1 ), 					label=L"y_r", linestyle=:dot, linewidth=2) 
		
	#Pinto los estados, si no quiero etiqueta
	fig2=plot(sol, xaxis="t",yaxis=L"x(t)", label=:none)
	#fig2=plot(sol,vars=(0,1), xaxis=L"t",yaxis=L"x(t)",label=L"x_1")
	#fig2=plot!(sol,vars=(0,2), ,label=L"x_2")
	    
	#Señal de control, u no se guarda, hay que recalcularlo
	fig3=plot(sol, vars=((t,x1,x2)->(t,control_1([x1,x2],parametros,t)),0, 1, 2),
		    xaxis=L"t",yaxis=L"u(t)",legend=false)
	l = @layout [a ; b ; c]
	figura=plot(fig1,fig2,fig3, layout=l)
	#display(figura) si se usa fuera de pluto hay que hacer display
end

# ╔═╡ 922bccde-022d-11eb-17f2-270aae61f923
md"""# Dinámica cero
El ejemplo anterior funciona muy bien por dos motivos:
1) Al derivar obtengo tantos estados como tenía antes (grado relativo máximo)
2) El cambio de variabes de los estados a las derivadas de y era invertible (difeomorfismo)
Eso significa que controlar $y$ y sus derivadas $y, \dot y ... y^{(\gamma)}$ es lo mismo que controlar $x_1...x_\gamma$.

Pero ¿y si no es así? ¿Qué pasa si u aparece "demasiado pronto"?

"""

# ╔═╡ cc0e02e2-04c9-11eb-2e08-b906aaec3959
md"""# Otro ejemplo

Segundo ejemplo

$$\dot x_1=-ax_1 + e^{2x_2}u$$
$$\dot x_2=2x_1x_2 +sin(x_2) + \frac{u}{2} $$
$$\dot x_3=2x_2$$
$$y=x_3$$
"""

# ╔═╡ 8839fb60-095a-11eb-063e-85fd960ea83a
md"Ya sabemos hacerlo, vamos al grano"

# ╔═╡ 54523da0-0958-11eb-2730-cbef0823f684
begin
	a=x3=0
	@variables a x3
	#copia de las ecuaciones poniendo Bs para no pisar variables anteriores
	dxb1=-x1 + exp(2x2)*u
    dxb2=2x1*x2 +sin(x2) + u/2
	dxb3=2x2	
	dxb=[dxb1 dxb2 dxb3]
	#estados y salida
	xb=[x1;x2;x3]
	yb=x3
end

# ╔═╡ 0610d2ee-0958-11eb-2490-19e0efcd3996
dyb=derivada(yb,xb,dxb,1)

# ╔═╡ 52eade10-095b-11eb-0717-a9686481d875
ddyb=expand(derivada(yb,xb,dxb,2))

# ╔═╡ 5c70004e-095b-11eb-0068-89e38689515b
md"""El grado relativo es 2, pero hay tres estados... *ups*

El jacobiano es obvio que es singular (2x3)
"""

# ╔═╡ ab5f1de0-0960-11eb-2cf6-414898fd1323
jacobian([yb;dyb],xb)

# ╔═╡ c7d09f80-0a28-11eb-3548-11dbb05bbbbb
md"""Vamos a seguir *como si nada* a ver que pasa...
"""

# ╔═╡ 77f78c80-a997-11eb-0b80-8da3a424a52f
string(solve_for(ddyb~v,u))

# ╔═╡ 706b0500-0a2d-11eb-023e-1de0809fe6cf
referencia_2(t)=0*referencia_1_a(t)

# ╔═╡ c8abde00-095c-11eb-237c-a9fd2bd95b8a
function control_2(x,p,t)
	x1,x2,x3=x
	
	yr,dyr,ddyr=referencia_2(t) #defindia más abajo
	
	d0=2
	d1=3
	
	y=x3
	dy=2*x2
	
	e0=yr-y
	e1=dyr-dy
	
	v=d0*e0 + d1*e1 + ddyr
	
	u= v - 2*sin(x2) - 4*x1*x2
	return u
end	

# ╔═╡ 97e9d100-04ca-11eb-07ea-3f5f5addbc45
function derivadas_ejemplo2(x,p,t)
	x1,x2,x3=x
	a=p[1]
	u=control_2(x,p,t);#calculamos el control
	dx1=-a*x1 + exp(2x2)*u
    dx2=2x1*x2 +sin(x2) + u/2
	dx3=2x2
	dx=[dx1 dx2 dx3]
end

# ╔═╡ c8a04540-095c-11eb-0b16-0bbf6ae0b6eb
md"""# ¡Funciona!
Pero si cambiamos a=1 por a=-1 ¿Qué pasaría?"""

# ╔═╡ 062714d0-0a2e-11eb-073f-6bd47b020c54
	constante_a=1

# ╔═╡ c900ad60-095f-11eb-3d18-9f187abdb9b3


# ╔═╡ 667f8b70-095f-11eb-08a8-0f49a920e3cb
@bind tf Slider(10:100)

# ╔═╡ b117ab50-04cc-11eb-1580-fdafb3c3213b
let
	x0 = [0.1 0.1 0.1]
	parametros=[constante_a]
	tspan = (0.0,tf)
	prob = ODEProblem(derivadas_ejemplo2,x0,tspan,parametros);
	sol = solve(prob,Tsit5());
	
	# Salida y referencia idxs=(0,3) significa pintar x3 frente a t
	fig1=plot(sol,idxs=(0,3), xaxis=L"t",yaxis=L"y(t)",label=L"y(t)")
	fig1=plot!(fig1,sol,idxs=( (t,x1)->(t,referencia_2(t)[1]) , 0, 1 ), 					label=L"y_r", linestyle=:dot, linewidth=2) 
		
	fig2=plot(sol, xaxis="t",yaxis=L"x(t)", label=:none)
	#fig2=plot(sol,idxs=(0,1), xaxis=L"t",yaxis=L"x(t)",label=L"x_1")
	#fig2=plot!(sol,idxs=(0,2), ,label=L"x_2")
	

	#Señal de control, u no se guarda, hay que recalcularlo
	fig3=plot(sol, idxs=((t,x1,x2, x3)->(t,control_2([x1,x2,x3],parametros,t)),0, 1, 2, 3),
		    xaxis=L"t",yaxis=L"u(t)",legend=false)
	l = @layout [a ; b ; c]
	figura=plot(fig1,fig2,fig3, layout=l)
	#display(figura) si se usa fuera de pluto hay que hacer display
end

# ╔═╡ 1f8be93e-0962-11eb-14f4-73a59ba80a72
md""" Que raro... para ver que ha pasado veamos la dinámica cero, recordemos las ecuaciones:

$\dot x_1=-ax_1 + e^{2x_2}u$
$\dot x_2=2x_1x_2 +sin(x_2) + \frac{u}{2}$
$\dot x_3=2x_2$
$y=x_3$

Si la salida $y(t)=0$ cero entonces:

$y=x_3=0 \to \dot x_3=0=2x_2 \to x_2=0 \to \dot x_2=0$
$0=0 +sin(0) + \frac{u}{2} \to u=0$

$\dot x_1=-ax_1$

Esto es estable si a es positivo pero *inestable* si a es negativo, eso podría explica por qué x_1 se aleja del origen (dinámica cero inestable)
"""

# ╔═╡ 41656a80-0a7c-11eb-1070-db649ce73728
md"""Para comprobar la teoría ponemos la referenica a cero, cuando la salida se haga cero lo que queda es "la dimámica cero" y la teoría nos dice:

-Dinámica cero asintóticamete estable $\to$ el sistema es **localmente** estable

-Dinámica cero inestable $\to$ el sistema es inestable

**OJO** lo que la teoría no dice es que pueda seguir una trayectoria, en este caso la sigue pero podría ser que no (al seguir una referencia *general* la salida no es cero y si es muy distinta de cero nos alejamos de la zona conocida)
"""

# ╔═╡ Cell order:
# ╟─5c4af6e0-022e-11eb-152f-2153b9b672c1
# ╠═85c0d350-022e-11eb-0873-c118ffde014c
# ╠═506a20d0-026f-11eb-2bde-f94724bd27f3
# ╟─85c03440-0274-11eb-3c81-514259334750
# ╟─2a98e6e0-a227-11eb-026a-bb5831f344bb
# ╟─e3e9d0a0-022c-11eb-04ba-97dcb9ba908a
# ╟─6831eb02-023b-11eb-0842-fbe8f8eb8937
# ╠═a2fa5231-07fe-448c-bbcb-a42bf790e7bc
# ╟─1c0442a0-a993-11eb-1b2d-55eb8b868dc7
# ╠═3ff84c30-0248-11eb-33ba-f99304ee8672
# ╠═6d535e50-a99a-11eb-08c1-4d98301e53b3
# ╠═5525d1a0-a9c7-11eb-1de0-1700c9c20460
# ╠═65893910-a9c7-11eb-2f9f-ebadf375cdbc
# ╠═a0571620-a9c7-11eb-070f-297fa9e406b2
# ╟─39ce8af0-023d-11eb-3178-857433b369ae
# ╟─4c478bf0-023d-11eb-0dce-439d408010fa
# ╟─b7e08f60-023d-11eb-13cb-ad5285e900d8
# ╠═697f24d0-a99a-11eb-26e6-a749fb04cb0f
# ╟─9ecc5ff0-022d-11eb-3715-51831818d04f
# ╠═65ad8360-a99a-11eb-3764-aded1015107e
# ╠═4b753d20-a9c8-11eb-320a-d58ec5e3c431
# ╠═6e31c180-a9c8-11eb-3893-f55415fa4fad
# ╠═478bc500-2c06-4508-8794-0d99fee916df
# ╟─561e07e2-0250-11eb-1380-4b5dfa861cd6
# ╠═64db3f40-a99a-11eb-01df-d3ca4846ca33
# ╟─de497f50-0250-11eb-028e-ffe56cf1ba45
# ╠═c2f97c40-a99a-11eb-30e7-359effec0a1d
# ╠═85e1f1f0-a9c9-11eb-2a97-1944f7fcb34b
# ╟─0db8c9c0-0252-11eb-199c-2d69997c6769
# ╠═dd56c0c2-a99a-11eb-0c87-8f8e2c6321a8
# ╟─72e4e740-022d-11eb-3a63-31023b9d464f
# ╠═f1e0ee80-a99a-11eb-0c06-85f8d7086aac
# ╠═03abf0e0-a9ca-11eb-339c-5599c951871c
# ╠═183a6460-a9ca-11eb-340e-872e08f29a9e
# ╠═2d9109e0-a9ca-11eb-08a5-7596eee9e4b2
# ╟─9ec274e0-022d-11eb-085a-7907b8f1eb93
# ╠═aed8d350-0265-11eb-1223-1dd0651632b6
# ╟─0d736c40-0266-11eb-3433-ad75c1617314
# ╠═454272d0-a9cb-11eb-36db-ffb885a95a52
# ╟─f54ff5b2-0266-11eb-1df3-49f3c95729c7
# ╠═2814f040-0267-11eb-3f4d-2933c0c3528a
# ╟─3b125190-0269-11eb-2473-33cbb7451779
# ╠═639dce00-0269-11eb-191e-d16a91f3bcf9
# ╟─e0b21ca0-0a2c-11eb-0e50-efa6c355e7c5
# ╟─f162f3f0-022d-11eb-30d9-55fb6e39c8c0
# ╠═a7745f70-026b-11eb-1abe-d1944165d231
# ╠═a0e6f960-026b-11eb-23c9-410e6d61c4a1
# ╟─922bccde-022d-11eb-17f2-270aae61f923
# ╟─cc0e02e2-04c9-11eb-2e08-b906aaec3959
# ╠═97e9d100-04ca-11eb-07ea-3f5f5addbc45
# ╟─8839fb60-095a-11eb-063e-85fd960ea83a
# ╠═54523da0-0958-11eb-2730-cbef0823f684
# ╠═0610d2ee-0958-11eb-2490-19e0efcd3996
# ╠═52eade10-095b-11eb-0717-a9686481d875
# ╟─5c70004e-095b-11eb-0068-89e38689515b
# ╠═ab5f1de0-0960-11eb-2cf6-414898fd1323
# ╟─c7d09f80-0a28-11eb-3548-11dbb05bbbbb
# ╠═77f78c80-a997-11eb-0b80-8da3a424a52f
# ╠═c8abde00-095c-11eb-237c-a9fd2bd95b8a
# ╠═706b0500-0a2d-11eb-023e-1de0809fe6cf
# ╟─c8a04540-095c-11eb-0b16-0bbf6ae0b6eb
# ╠═062714d0-0a2e-11eb-073f-6bd47b020c54
# ╠═c900ad60-095f-11eb-3d18-9f187abdb9b3
# ╠═667f8b70-095f-11eb-08a8-0f49a920e3cb
# ╠═b117ab50-04cc-11eb-1580-fdafb3c3213b
# ╟─1f8be93e-0962-11eb-14f4-73a59ba80a72
# ╟─41656a80-0a7c-11eb-1070-db649ce73728
