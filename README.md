# Guia-de-instalação-do-FEniCSx-para-Windows
Um guia para a instalação do [FEniCSx 0.10v](https://fenicsproject.org/) para Windows 11. Nesse guia estarei assumindo que você já possua o [python3](https://www.python.org/downloads/) e algum ambiente de desenvolvimento, como o [VSCode](https://code.visualstudio.com/docs/?dv=win64user), instalados. 

# Windows Subsystem for Linux (WSL)

Primeiramente, o FEniCSx é uma plataforma que possui suporte nativo para certas distribuições Linux. Uma delas, e a que iremos utilizar, será o Ubuntu 24.04+.
Por sorte, o Windows possui um recurso que permite rodar um ambiente Linux dentro dele sem precisar realizar dual-boot ou uma VM. Esse recurso é o [WSL](https://learn.microsoft.com/pt-br/windows/wsl/install).

Para realizar sua instalação, será necessário abrir o PowerShell no modo administrador.

```
wsl --install
```

Após habilitado o WSL, iremos baixar o Ubuntu
```
wsl --install Ubuntu-24.04
```

Ao copiar o código acima, você realizará a instalação da distro. Será necessário definir o nome e a senha para o usuário sudo do novo sistema, então escolha uma senha de fácil acesso na qual você se lembrará depois.

![image](https://github.com/SSputnik82/Guia-de-instala-o-do-FEniCSx-para-Windows/blob/main/Captura%20de%20tela%202026-01-05%20103628.png)

>[!NOTE]
>Tenha certeza que você está usando a versão WSL2. Para saber qual versão você está utilizando copie o seguinte comando no PowerShell.
> ```
> wsl --status
> ```
> Verifique se a versão padrão é 2.
>![Image](https://github.com/SSputnik82/Guia-de-instala-o-do-FEniCSx-para-Windows/blob/main/Captura%20de%20tela%202026-01-05%20105445.png)

# Instalação do FEniCSx

Agora, abra o terminal do Ubuntu, você pode pesquisar por ele na área de pesquisa do próprio Windows, ou pode copiar o seguinte comando no PowerShell.

```
wsl ~ -d Ubuntu
```

Após entrar no ambiente Linux, iremos iniciar a instalação do FEniCSx nele. Copie os seguintes comandos no terminal em ordem.

```
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt update
sudo apt install fenicsx
```

conforme você for colando as linhas de comando no terminal, ele irá pedir a senha sudo criada anteriormente. Basta digitar ela e apertar Enter. As vezes irá aparecer a opção de [Y/N] que significa sim e não, basta digitar "Y" e seguir com a instalação.

# Conexão com o VSCode

Agora precisamos 
