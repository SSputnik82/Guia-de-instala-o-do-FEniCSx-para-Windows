# Guia de instalação do FEniCSx para Windows
Um guia para a instalação do [FEniCSx 0.10v](https://fenicsproject.org/) para Windows 11. Nesse guia estarei assumindo que você já possua o [python3](https://www.python.org/downloads/) e algum ambiente de desenvolvimento, como o [VSCode](https://code.visualstudio.com/docs/?dv=win64user), instalados. 

# Windows Subsystem for Linux (WSL)

Primeiramente, o FEniCSx é uma plataforma que possui suporte nativo para certas distribuições Linux. Uma delas, e a que iremos utilizar, será o Ubuntu 24.04+.
Por sorte, o Windows possui um recurso que permite rodar um ambiente Linux dentro dele sem precisar realizar dual-boot ou uma VM. Esse recurso é o [WSL](https://learn.microsoft.com/pt-br/windows/wsl/install).

Para realizar sua instalação, será necessário abrir o PowerShell no modo administrador.

```
wsl --install
```

Após habilitado o WSL, reinicie o computador e em seguida baixe o Ubuntu utilizando o comando:
```
wsl --install Ubuntu-24.04
```

Ao copiar o código acima, você realizará a instalação da distro. Será necessário definir o nome e a senha para o usuário sudo do novo sistema, então escolha uma senha de fácil acesso na qual você se lembrará depois.
>[!NOTE]
>Importante destacar que no momento em que você colocar a senha, ela aparecerá invisível. Não se preocupe e continue digitando normalmente pois ela está sendo escrita.
>Caso você cometa um erro e queria refazer a instalação, copie e cole o comando no PowerShell:
>```
>wsl --unregister Ubuntu-24.04
>```

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

Agora você poderá rodar seus códigos digitando no terminal da seguinte forma:
```
python3 nome_do_seu_codigo.py
```

<img width="1546" height="213" alt="Captura de tela 2026-01-05 140911" src="https://github.com/user-attachments/assets/d77dfa12-6bb3-4662-b6c7-25d2674d1241" />

## Troubleshooting
>[!TIP]
>Um dos possíveis problemas que podem aparecer é:
><img width="1919" height="161" alt="Captura de tela 2026-01-05 212142" src="https://github.com/user-attachments/assets/bc18aa75-36c3-402b-af2e-f17dd4fd9fe5" />
>
>Esse erro significa que seu terminal não está no mesmo diretório, ou pasta, que o seu arquivo. Para corrigir isso é preciso que você procure onde ele foi salvo no Windows. 
>```
>C:\Users\tsuno\Desktop\FEniCSx
>```
>No Ubuntu, utilizamos o comando <mark> cd </mark> para nos movermos entre diretórios. Porém, como nosso arquivo existe no sistema Windows, a forma de conectar ambos os ambientes é utilizando o comando <mark> mnt </mark>. Estaremos "montando" o diretório do Windows no sistema Linux. Lembrando que o diretório muda de usuário para usuário, portanto você deverá digitar no seu terminal o diretório pertinente à sua pasta:
>```
>cd /mnt/c/Users/tsuno/Desktop/FEniCSx
>```
><img width="1917" height="191" alt="Captura de tela 2026-01-05 213746" src="https://github.com/user-attachments/assets/8e8b37f0-7b58-4357-a398-d825fda333d0" />



# ParaView
Ele é um aplicativo de código aberto para visualização interativa das simulações que realizaremos no FEniCSx. Sua maior vantajem é a sua capacidade de permitir a computação paralela para problemas que requerem alta demanda computacional. Para realizar sua instalação, basta ir em seu site, ir para Downloads e baixar a versão no formato <mark> AMD64.msi </mark> mais recente: [ParaView](https://www.paraview.org/).
<img width="1919" height="657" alt="Captura de tela 2026-01-05 210111" src="https://github.com/user-attachments/assets/23c1885a-6682-493e-b61a-3a201dc2369b" />



Para abrir um arquivo a ser visualisado, basta apertar <mark> ctrl + o </mark> e procurar o diretório do seu arquivo.

<img width="1919" height="994" alt="Captura de tela 2026-01-05 142347" src="https://github.com/user-attachments/assets/0e408905-95c7-4303-ba77-87a732e95e79" />

> [!NOTE]
> O paraview consegue reconhecer dados por meio de uma variedade de formatos de arquivos. Porém, iremos focar apenas em alguns para fins de praticidade.
> Eles são <mark> .xdmf </mark>, <mark> .pvd </mark> e <mark> .vtm </mark>. Sendo esses dois últimos parte do formato VTK.

# Código teste
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
E = 1e+4
c = 1.0
A = 1.0
Le = 1.0

#Malha
domain = create_interval(MPI.COMM_WORLD, nx=32, points=(0,2*Le))
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
t_bar = Constant(domain, default_scalar_type(-c*Le**2 / A))

#a = E*u.dx(0)*v.dx(0)*dx
a = ufl.dot(E*ufl.grad(u),ufl.grad(v))*ufl.dx
L = b*v*dx + t_bar*v*ds

problem = LinearProblem(a, L, bcs=[bc_D],
    petsc_options_prefix="basic_linear_problem",
    petsc_options= {
      "ksp_type": "preonly",
      "pc_type": "lu",
      "pc_factor_mat_solver_type": "mumps"
})
uh = problem.solve()

xdmf = io.XDMFFile(domain.comm, "bar_el_1D.xdmf", "w")
xdmf.write_mesh(domain)
xdmf.write_function(uh)
xdmf.close()


```
Ao rodar o código seu terminal ficará assim:
<img width="1919" height="140" alt="Captura de tela 2026-01-05 210351" src="https://github.com/user-attachments/assets/298b5c59-df22-4cab-8e39-a9a145f7ff02" />

Agora abra o ParaView, aperte <mark> ctrl + o </mark> e procure o arquivo salvo com nome de <mark> bar_el_1D.xdmf </mark>. Em seguida escolha um dos leitores de arquivo e clique em <mark> Apply </mark>


<img width="1919" height="989" alt="Captura de tela 2026-01-05 192900" src="https://github.com/user-attachments/assets/f01e96e1-6653-49aa-9f1e-70571a28ba77" />


