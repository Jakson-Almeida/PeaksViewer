# PeaksViewer

Visualizador interativo de espectros ópticos com detecção de picos e ajuste de curvas, em Python (Tkinter + Matplotlib).

Permite carregar um ou vários arquivos de espectro (formato `wavelength;intensity`), navegar entre eles com o teclado e analisar picos, larguras (FWHM) e perfis Gaussianos/Lorentzianos.

## Funcionalidades

- **Carregamento múltiplo**: abra vários arquivos de uma vez e navegue entre eles.
- **Navegação por teclado**: `←` / `→`, `<` / `>`, `,` / `.`, `-` / `+`, `Page Up` / `Page Down`.
- **Detecção automática de picos** via `scipy.signal.find_peaks`, com prominência ajustável (filtra ruído).
- **Clique em pico** para selecioná-lo e copiar λ (comprimento de onda) ou I (intensidade) para a área de transferência.
- **Gradiente de cores**: preenchimento sob a curva colorido segundo o comprimento de onda (espectro visível 380–720 nm).
- **Escala em dB**: alterna entre intensidade linear e potência em dB (relativa ao máximo).
- **Ajuste de curva** Gaussiano ou Lorentziano com cálculo de:
  - Centro (λ₀)
  - Amplitude
  - FWHM (largura à meia altura)
  - R² (qualidade do ajuste)
- **Seleção de região** com o mouse (arrastar sobre o gráfico) para ajustar apenas um intervalo do espectro.

## Formato dos arquivos de entrada

Arquivos texto (`.txt` ou `.csv`), uma amostra por linha, separados por `;`:

```
6.5000000e-07;120.5
6.5010000e-07;125.2
6.5020000e-07;131.8
...
```

- **Coluna 1**: comprimento de onda em **metros** (será convertido para nm internamente: `× 1e9`).
- **Coluna 2**: intensidade (unidade arbitrária).

Linhas inválidas são ignoradas silenciosamente.

## Instalação

Requer Python 3.8+ e as seguintes bibliotecas:

```bash
pip install numpy scipy matplotlib
```

> `tkinter` faz parte da biblioteca padrão do Python (no Windows e macOS já vem incluído; no Linux pode ser necessário instalar `python3-tk`).

## Uso

```bash
python peaks_viewer.py
```

1. Clique em **Carregar arquivo(s)...** e selecione um ou mais espectros.
2. Use os botões `< Anterior` / `Próximo >` ou as teclas de seta para navegar.
3. Marque **Exibir picos** para detectar picos automaticamente; ajuste a **Prominência** para filtrar ruído.
4. Marque **Gradiente de cores** para preencher a área sob a curva com as cores do espectro visível.
5. Marque **Potência (dB)** para visualizar em escala logarítmica.
6. Marque **Ajustar curva**, escolha **gaussian** ou **lorentzian** e, opcionalmente, arraste o mouse sobre o gráfico para restringir o ajuste a uma região.
7. Clique em um pico e use **Copiar λ** ou **Copiar I** para copiar os valores.

## Atalhos de teclado

| Tecla                  | Ação              |
| ---------------------- | ----------------- |
| `←`, `<`, `,`, `-`, `Page Up`   | Espectro anterior |
| `→`, `>`, `.`, `+`, `Page Down` | Próximo espectro  |

## Estrutura

- `peaks_viewer.py` — aplicação completa (GUI, leitura de arquivos, detecção de picos, ajuste de curvas).

## Licença

Sem licença especificada.
