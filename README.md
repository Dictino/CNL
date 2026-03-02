# CNL
Herramientas del curso de [CNL](http://portal.uned.es/portal/page?_pageid=93,70656202&_dad=portal&_schema=PORTAL&idAsignatura=31104178&idTitulacion=310401) del [Máster en Ingeniería de Sistemas y Control UNED/UCM](https://cv4.ucm.es/moodle/course/view.php?id=4056)

## Instalación
Para ejecutar los ejemplos lo primero que necesitamos es instalar Julia, para ello lo mejor es ir a https://julialang.org/ y visitar su página de [descargas](https://julialang.org/downloads/) y seguir las instrucciones para instalar la última versión estable.

A continuación, desgargad los ficheros de este repositorio en vuestro ordenador y hecho esto desde un terminal:

```
$ julia "C:\ruta\en\mi\ordenador\CNL.jl"

```

También podéis hacerlo desde dentro de julia:

```julia
julia> include("C://ruta//en//mi//ordenador//CNL.jl")

```


*NOTA:* podéis usar la barra espaciadora para ir autocompletando la dirección o copiar y pegar desde el explorador y cambiar el separdor de directorios por //. También es cómodo pulsar la flecha de arriba o abajo para navegafr por la historia de órdenes y no tener que estar buscando la ruta o incluso crear un acceso directo con dicho comando.


Esto **tardará un buen rato** la primera vez ya que tiene que decargar e instalar todos los paquetes que usará el curso (Pluto,  Plots, OrdinaryDiffEq ModelingToolkit, SymEngine...) así que aprovechad para ver algún tema de CNL ;)

Si se quejase por un problema de versiones probad a borrar "Manifest.toml" e "instalado" y volved a intentarlo.

Cuando termine se abrirá pluto en el navegador:

![pluto](https://github.com/Dictino/CNL/blob/main/Im%C3%A1genes/pluto.png?raw=true)

Desde ahí se puede seleccionar el notebook que se desea abrir pinchando sobre la barra, buscando el fichero y pulsando Open.

*NOTA:* No es necesario usar Pluto ni este método para leer o utilizar el código, los notebooks .jl son scripts normales de julia y podrían ejecutarse de la misma forma que CNL.jl y abrise con cualquier editor de textos. Eso sí, se perdería la interactividad, y para mostrar los resultados había cambiar el código usar ```@show variable``` o ```println("texto que sea"")``` y para los graficos habría que ```display(figura_correspondiente)```.

¡A disfrutarlo!

## Primeros pasos con Julia

En los vídeos de la asignatura os explciaremos todo lo que necesitais conocer de la sintaxis del lenguaje pero para ampliar conocimientos o resolver dudas es bueno contar con los siguietnes recursos:

En primer lugar una guía para conocer la sintaxis básica de Julia es la siguiente [tutorial en castellano](https://hedero.webs.upv.es/julia-basico/), o también [estos en inglés](https://julialang.org/learning/tutorials/) (el primero "From zero to Julia!" merece bastante la pena).

Un recurso fundamental para la consulta y también para leer secciones completas ya que es bastante legible para ser un manual es la [documentación oficial](https://docs.julialang.org/en/v1/) que es muy completa y exausiva y si se quiere concoer el lenguaje en profundidad el manual es de obligada lectura.

Para los que dominéis otros lenguajes puede ser interesante la siguiente [sección](https://docs.julialang.org/en/v1/manual/noteworthy-differences/) donde explica las diferencias con otros lenguajes de amplio uso.

Finalmente cada paquete tiene su propia documentación que puede consultarse buscando en google "julia paquete_que_sea.jl" por ejemplo "julia differentialequations.jl" en el primer resultado lleva directamente a la [documentación](https://docs.sciml.ai/DiffEqDocs/stable/) de dicho paquete.
 
