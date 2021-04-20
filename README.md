# CNL
Herramientas del curso de [CNL](http://portal.uned.es/portal/page?_pageid=93,70656202&_dad=portal&_schema=PORTAL&idAsignatura=31104178&idTitulacion=310401) del [Máster en Ingeniería de Sistemas y Control UNED/UCM](https://cv4.ucm.es/moodle/course/view.php?id=4056)

## Instalación
Para ejecutar los ejemplos lo primero que necesitamos es instalar Julia, para ello lo mejor es ir a https://julialang.org/ y visitar su página de [descargas](https://julialang.org/downloads/) y descargar a última versión estable. La instalación es trivial, basta con ejecutar el fichero descargado o en el caso de Linux basta con descomprimir el fichero tar en una carpeta cualquiera de vuestra elección y después ejecutar julia que se encuentra en la subcarpeta /bin/julia. (Para comodidad es conveniente crear un enlace simbólico desde la carpeta donde se encuentra julia a una carpeta que esté en el path ej. ```sudo ln -s $(pwd)/julia /bin/julia```)

A continuación, desde el terminal de julia:

```julia
julia>] add Pluto PlutoUI Plots OrdinaryDiffEq ModelingToolkit SymEngine LaTeXStrings ControlSystems Symbolics

```

Esto tardará un rato (sobre todo con mala conexión a internet y si se usa windows) así que aprovechad para ver algún tema de CNL.

Cuando acabe podéis lanzar pluto:

```julia
julia> using Pluto
julia> Pluto.run()
```

Y se abrirá pluto en el navegador:

![pluto](https://github.com/Dictino/CNL/blob/main/Im%C3%A1genes/pluto.png?raw=true)

Desde ahí se puede seleccionar el notebook que se desea abrir, o bien descargando los ficheros en vuesto ordenador o poniendo por ejemplo la url de github de dicho notebook:

https://github.com/Dictino/CNL/blob/main/Tema%201.jl

https://github.com/Dictino/CNL/blob/main/Tema%202.jl

https://github.com/Dictino/CNL/blob/main/Tema%203.jl

https://github.com/Dictino/CNL/blob/main/Tema%204.jl

https://github.com/Dictino/CNL/blob/main/Tema%205.jl

https://github.com/Dictino/CNL/blob/main/Tema%206.jl

Si queréis probarlo **sin instalar nada** en vuestro ordenador (no es la opción recomendada, pero es muy útil si queréis verlo en una tableta, teléfono o en un ordenador en el que no queréis instalar software)


[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Dictino/pluto-on-binder/2a61735?urlpath=pluto)

Ojo que binder **no es persistente**, si hacéis un cambio y queréis guardarlo no uséis el botón de guardar sino el de exportar:

![pluto](https://github.com/Dictino/CNL/blob/main/Im%C3%A1genes/Boton_exportar.png?raw=true)

¡A disfrutarlo!

Una guía para conocer la sintaxis básica de Julia es el siguiente [tutorial en castellano](https://hedero.webs.upv.es/julia-basico/), o también [estos en inglés](https://julialang.org/learning/tutorials/) (el primero "From zero to Julia!" merece bastante la pena).

Un recurso fundamental es la [documentación](https://docs.julialang.org/en/v1/) que es muy completa y exausiva y si se quiere concoer el lenguaje el manual es de obligada lectura, para los que dominéis otros lenguajes puede ser interesante la [sección](https://docs.julialang.org/en/v1/manual/noteworthy-differences/).
 
