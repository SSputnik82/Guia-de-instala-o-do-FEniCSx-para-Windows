# Guia de instalação do FEniCSx para Windows
Um guia para a instalação do [FEniCSx 0.10v](https://fenicsproject.org/) para Windows 11. Nesse guia estarei assumindo que você já possua o [python3](https://www.python.org/downloads/) e algum ambiente de desenvolvimento, como o [VSCode](https://code.visualstudio.com/docs/?dv=win64user), instalados. 

# Windows Subsystem for Linux (WSL)

Primeiramente, o FEniCSx é uma plataforma que possui suporte nativo para certas distribuições Linux. Uma delas, e a que iremos utilizar, será o Ubuntu 24.04+.
Por sorte, o Windows possui um recurso que permite rodar um ambiente Linux dentro dele sem precisar realizar dual-boot ou uma VM. Esse recurso é o [WSL](https://learn.microsoft.com/pt-br/windows/wsl/install).

Para realizar sua instalação, será necessário abrir o PowerShell no modo administrador.

```
wsl --install
```

Após habilitado o WSL, reinicie o PowerShell e em seguida baixe o Ubuntu utilizando o comando:
```
wsl --install Ubuntu-24.04
```

Ao copiar o código acima, você realizará a instalação da distro. Será necessário definir o nome e a senha para o usuário sudo do novo sistema, então escolha uma senha de fácil acesso na qual você se lembrará depois.

![image](https://github.com/SSputnik82/Guia-de-instala-o-do-FEniCSx-para-Windows/blob/main/Captura%20de%20tela%202026-01-05%20103628.png)

>[!NOTE]
>Tenha certeza que você está usando a versão WSL2. Para saber qual a sua versão, copie o seguinte comando no PowerShell.
> ```
> wsl --status
> ```
> Verifique se a versão padrão é a 2, caso não seja defina como padrão:
>```
>wsl --set-default-version 2
>```
>![Image](https://github.com/SSputnik82/Guia-de-instala-o-do-FEniCSx-para-Windows/blob/main/Captura%20de%20tela%202026-01-05%20105445.png)

# Instalação do FEniCSx

Agora, abra o terminal do Ubuntu, você pode pesquisar por ele na área de pesquisa do próprio Windows, ou pode copiar o seguinte comando no PowerShell.

```
wsl ~ -d Ubuntu
```

Após entrar no ambiente Linux, iremos iniciar a instalação do FEniCSx nele. Copie e cole os seguintes comandos, um de cada vez, no terminal.

```
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt update
sudo apt install fenicsx
```

Conforme você for colando as linhas de comando no terminal, ele irá pedir a senha sudo criada anteriormente. Basta digitar ela e apertar Enter. Às vezes irá aparecer a opção de [Y/N] que significa sim e não, basta digitar "Y" e apertar Enter para seguir com a instalação. Após terminada reinicie o computador.

# Integração com o VSCode

Para utilizar o terminal do Ubuntu no VSCode, basta abrir <mark> View > Terminal </mark>, ou clicar no botão de <mark> Run python file </mark>.
<img width="1919" height="1017" alt="Captura de tela 2026-01-05 141144" src="https://github.com/user-attachments/assets/4a81b072-2adc-43bf-9483-e0d2af08ff4b" />

Após abrir o terminal, procure a opção <mark> Ubuntu (WSL) </mark>.
<img width="1919" height="1021" alt="Captura de tela 2026-01-05 144825" src="https://github.com/user-attachments/assets/b2a375e6-1693-419d-beb0-ec45c3b41f29" />

<img width="1546" height="213" alt="Captura de tela 2026-01-05 140911" src="https://github.com/user-attachments/assets/d77dfa12-6bb3-4662-b6c7-25d2674d1241" />

Agora você poderá rodar seus códigos digitando no terminal o comando.
```
python3 Nome_do_codigo.py
```


# ParaView
Ele é um aplicativo de código aberto para visualização interativa das simulações que realizaremos no FEniCSx. Sua maior vantajem é a sua capacidade de permitir a computação paralela para problemas que requerem alta demanda computacional. Para realizar sua instalação, basta ir em seu site, ir para Downloads e baixar a versão mais recente: [ParaView](https://www.paraview.org/).

Para abrir um arquivo a ser visualisado, basta apertar <mark> ctrl + o </mark> e procurar o diretório do seu arquivo.

<img width="1919" height="994" alt="Captura de tela 2026-01-05 142347" src="https://github.com/user-attachments/assets/0e408905-95c7-4303-ba77-87a732e95e79" />

> [!NOTE]
> O paraview consegue reconhecer dados por meio de uma variedade de formatos de arquivos. Porém, iremos focar apenas em alguns desses formatos para fins de praticidade.
> Eles são <mark> .xdmf </mark>, <mark> .pvd </mark> e <mark> .vtm </mark>. Sendo esses dois últimos parte do formato VTK.

Código teste:
```python

#Barra elastica 1D sob tração/compressão

import ufl

from dolfinx import default_scalar_type, io
from dolfinx.fem import (Constant, Function, functionspace,
						 dirichletbc, locate_dofs_geometrical)
from dolfinx.fem.petsc import LinearProblem
from dolfinx.mesh import create_interval

from mpi4py import MPI
from ufl import SpatialCoordinate, TestFunction, TrialFunction, dot, dx, ds, grad

import numpy as np

#Variaveis
E = 1E+4
c = 1.0
A = 1.0
L= 1.0

#Malha
domain = create_interval(MPI.COMM_WORLD, nx=32, points=(0,2*L))
V = functionspace(domain, ("Lagrange", 1))
u = TrialFunction(V)
v = TestFunction(V)


#Condições de contorno
def boundary_D(x):
	return np.isclose(x[0],0)

dofs_D = locate_dofs_geometrical(V,boundary_D)
bc_D = dirichletbc(default_scalar_type(0), dofs_D, V)

#dx = ufl.Measure("dx", domain=domain)
x = SpatialCoordinate(domain)
b = c*x[0]
t_bar = Constant(domain, default_scalar_type(-c*L**2 / A))

#a = E*u.dx(0)*v.dx(0)*dx
a = ufl.dot(E*ufl.grad(u),ufl.grad(v))*ufl.dx
L = b*v*dx + t_bar*v*ds

problem = LinearProblem(a,L,bcs=[bc_D], petsc_options={"ksp_type": "preonly", "pc_type":"lu"})
uh = problem.solve()

xdmf = io.XDMFFile(domain.comm, "bar_el_1D.xdmf", "w")
xdmf.write_mesh(domain)
xdmf.write_function(uh)
xdmf.close()


```
Resultado visualisado no ParaView:
<img width="1919" height="989" alt="Captura de tela 2026-01-05 192900" src="https://github.com/user-attachments/assets/f01e96e1-6653-49aa-9f1e-70571a28ba77" />


