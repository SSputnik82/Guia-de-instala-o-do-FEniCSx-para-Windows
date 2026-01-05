# Guia de instalação do FEniCSx para Windows
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

conforme você for colando as linhas de comando no terminal, ele irá pedir a senha sudo criada anteriormente. Basta digitar ela e apertar Enter. Às vezes irá aparecer a opção de [Y/N] que significa sim e não, basta digitar "Y" e apertar Enter para seguir com a instalação.

# Integração com o VSCode

Para utilizar o terminal do Ubuntu no VSCode, basta abrir <mark> View > Terminal </mark>, ou clicar no botão de <mark> Run python file </mark>.
(<img width="1919" height="1017" alt="Captura de tela 2026-01-05 141144" src="https://github.com/user-attachments/assets/4a81b072-2adc-43bf-9483-e0d2af08ff4b" />)

Após abrir o terminal, procure a opção <mark> Ubuntu (WSL) </mark>.
<img width="1919" height="1022" alt="Captura de tela 2026-01-05 140524" src="https://github.com/user-attachments/assets/e803f5f5-be1b-486a-99ef-823e1d795ae9" />

Agora você poderá rodar seus códigos digitando no terminal o comando.
```
python3 [Nome_do_codigo.py]
```
<img width="1546" height="213" alt="Captura de tela 2026-01-05 140911" src="https://github.com/user-attachments/assets/1988fe82-b564-4cdf-9fe6-c65330cda1b8" />


# ParaView
Ele é um aplicativo de código aberto para visualização interativa das simulações que realizaremos no FEniCSx. Sua maior vantajem é a sua capacidade de permitir a computação paralela para problemas que requerem alta demanda computacional. Para realizar sua instalação, basta ir em seu site, ir para Downloads e baixar a versão mais recente: [ParaView](https://www.paraview.org/).

Para abrir um arquivo a ser visualisado, basta apertar <mark> ctrl + o </mark> e procurar o diretório do seu arquivo.

<img width="1919" height="994" alt="Captura de tela 2026-01-05 142347" src="https://github.com/user-attachments/assets/0e408905-95c7-4303-ba77-87a732e95e79" />

> [!NOTE]
> O paraview consegue reconhecer dados por meio de uma variedade de formatos de arquivos. Porém, iremos focar apenas em alguns desses formatos para fins de praticidade.
> Eles são <mark> .xdmf </mark>, <mark> .pvd </mark> e <mark> .vtm </mark>. Esses dois últimos sendo parte do formato VTK.

