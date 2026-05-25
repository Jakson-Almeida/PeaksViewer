"""
Visualizador de espectros com detecção de picos.
Permite carregar um ou vários arquivos de espectro (formato wl;intensity),
visualizar um por vez e navegar entre eles com < e > (ou setas).
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.collections import PolyCollection
from matplotlib.widgets import SpanSelector
import os


def wavelength_to_rgb(wavelength, gamma=0.8, dark=False):
    """
    Converte comprimento de onda (nm) para RGB.
    Espectro visível (420–680 nm): cores do arco-íris.
    Fora do intervalo: preto, com degradê entre a última cor e o preto.
    """
    wavelength = float(wavelength)
    # Abaixo de 380 nm: preto
    if wavelength < 380:
        return (0.0, 0.0, 0.0)
    # 380–420 nm: degradê de preto até a cor no limite (violeta em 420 nm)
    if wavelength < 420:
        t = (wavelength - 380) / (420 - 380)
        r = t ** gamma
        g = 0.0
        b = t ** gamma
        return (r, g, b)
    # Acima de 720 nm: preto
    if wavelength >= 720:
        return (0.0, 0.0, 0.0)
    # 680–720 nm: degradê da última cor (vermelho em 680 nm) até preto
    if wavelength > 680:
        t = (720 - wavelength) / (720 - 680)
        r = t ** gamma
        g = 0.0
        b = 0.0
        return (r, g, b)
    # Espectro visível 420–680 nm
    if 420 <= wavelength < 440:
        r = ((-(wavelength - 440) / (440 - 420))) ** gamma
        g = 0.0
        b = 1.0 ** gamma
    elif 440 <= wavelength < 490:
        r = 0.0
        g = ((wavelength - 440) / (490 - 440)) ** gamma
        b = 1.0 ** gamma
    elif 490 <= wavelength < 510:
        r = 0.0
        g = 1.0 ** gamma
        b = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif 510 <= wavelength < 580:
        r = ((wavelength - 510) / (580 - 510)) ** gamma
        g = 1.0 ** gamma
        b = 0.0
    elif 580 <= wavelength < 645:
        r = 1.0 ** gamma
        g = (-(wavelength - 645) / (645 - 580)) ** gamma
        b = 0.0
    else:  # 645 <= wavelength <= 680
        r = 1.0 ** gamma
        g = 0.0
        b = 0.0
    return (max(0, r), max(0, g), max(0, b))


def precompute_gradient(wl_nm, dark=False):
    """Pré-calcula as cores do gradiente para cada segmento do espectro (wl em nm)."""
    return [wavelength_to_rgb((wl_nm[j] + wl_nm[j + 1]) / 2, dark=dark) for j in range(len(wl_nm) - 1)]


def gaussian(x, amp, center, sigma):
    """Função gaussiana: amp * exp(-(x-center)^2 / (2*sigma^2))"""
    return amp * np.exp(-((x - center) ** 2) / (2 * sigma ** 2))


def lorentzian(x, amp, center, gamma):
    """Função lorentziana: amp * gamma^2 / ((x-center)^2 + gamma^2)"""
    return amp * gamma ** 2 / ((x - center) ** 2 + gamma ** 2)


def ajustar_curva(wl_nm, spec, modelo="gaussian"):
    """
    Ajusta uma curva gaussiana ou lorentziana aos dados do espectro.
    
    Args:
        wl_nm: Array de comprimentos de onda (nm)
        spec: Array de intensidades
        modelo: "gaussian" ou "lorentzian"
    
    Returns:
        (params, curva_ajustada, r_squared, fwhm) ou (None, None, None, None) se falhar
        params: (amp, center, width) onde width é sigma (gaussian) ou gamma (lorentzian)
    """
    try:
        # Estimativa inicial dos parâmetros
        amp_guess = np.max(spec)
        center_guess = wl_nm[np.argmax(spec)]
        width_guess = (wl_nm[-1] - wl_nm[0]) / 10  # ~10% da faixa
        
        p0 = [amp_guess, center_guess, width_guess]
        
        # Ajuste
        if modelo == "gaussian":
            popt, _ = curve_fit(gaussian, wl_nm, spec, p0=p0, maxfev=5000)
            curva = gaussian(wl_nm, *popt)
            fwhm = 2.355 * abs(popt[2])  # FWHM = 2.355 * sigma
        else:  # lorentzian
            popt, _ = curve_fit(lorentzian, wl_nm, spec, p0=p0, maxfev=5000)
            curva = lorentzian(wl_nm, *popt)
            fwhm = 2 * abs(popt[2])  # FWHM = 2 * gamma
        
        # Calcula R²
        ss_res = np.sum((spec - curva) ** 2)
        ss_tot = np.sum((spec - np.mean(spec)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        return popt, curva, r_squared, fwhm
    
    except Exception:
        return None, None, None, None


def ler_dados_arquivo(caminho_arquivo):
    """
    Lê arquivo de espectro (uma linha por ponto: wavelength;intensity).
    Retorna (wl_metros, intensidade) ou ([], []) em caso de erro.
    """
    frequencias = []
    ganhos = []
    try:
        with open(caminho_arquivo, "r", encoding="utf-8", errors="replace") as arquivo:
            for linha in arquivo:
                dados = linha.strip().split(";")
                if len(dados) >= 2:
                    try:
                        frequencias.append(float(dados[0]))
                        ganhos.append(float(dados[1]))
                    except ValueError:
                        continue
    except FileNotFoundError:
        return [], []
    except Exception:
        return [], []
    return frequencias, ganhos


def detectar_picos(intensidade, prominence=5, valley=False):
    """Retorna índices dos picos (ou vales) e (wl_nm, intensity) nesses pontos."""
    arr = np.asarray(intensidade)
    if valley:
        peaks, _ = find_peaks(-arr, prominence=prominence)
    else:
        peaks, _ = find_peaks(arr, prominence=prominence)
    return peaks


def _to_db(intensity, ref=None):
    """Converte intensidade para dB (ref = max se None). Evita log(0) com piso."""
    arr = np.asarray(intensity, dtype=float)
    if ref is None:
        ref = np.max(arr)
    if ref <= 0:
        ref = 1.0
    return 10 * np.log10(np.maximum(arr / ref, 1e-12))


def plotar_espectro_com_picos(ax, wl_nm, spec, prominence=5, valley=False, dark=False, show_peaks=False, show_gradient=False, fit_curve=None, fit_data=None, fit_wl=None, selected_range=None, power_db=False):
    """
    Plota espectro em ax. Se show_peaks=True, detecta picos e desenha marcadores.
    Se show_gradient=True, preenche a área sob a curva com gradiente de cores (λ).
    Se fit_curve e fit_data são fornecidos, plota a curva ajustada.
    Se selected_range é fornecido, desenha a região selecionada.
    Se power_db=True, eixo Y em dB (relativo ao máximo).
    Limpa marcadores antigos em ax (ax.markers e ax.marker).
    """
    ax.clear()
    color_fg = "white" if dark else "black"
    color_bg = "black" if dark else "white"
    ax.set_facecolor(color_bg)
    ax.set_xlabel("Comprimento de onda (nm)", color=color_fg)
    ax.set_ylabel("Potência (dB)" if power_db else "Intensidade (u.a.)", color=color_fg)
    ax.tick_params(colors=color_fg)
    ax.grid(True)

    ref = np.max(spec) if np.max(spec) > 0 else 1.0
    spec_plot = _to_db(spec, ref) if power_db else spec
    fit_data_plot = _to_db(fit_data, ref) if (power_db and fit_data is not None) else fit_data

    if show_gradient:
        gradient_colors = precompute_gradient(wl_nm, dark=dark)
        floor = np.min(spec_plot) - 10 if power_db else 0
        verts = [
            [(wl_nm[j], floor), (wl_nm[j], spec_plot[j]), (wl_nm[j + 1], spec_plot[j + 1]), (wl_nm[j + 1], floor)]
            for j in range(len(wl_nm) - 1)
        ]
        poly = PolyCollection(verts, facecolors=gradient_colors, edgecolors="none")
        ax.add_collection(poly)
        ax.plot(wl_nm, spec_plot, color="white" if dark else "gray", lw=1.5, alpha=0.8, label="Dados")
    else:
        ax.plot(wl_nm, spec_plot, color="gray" if dark else "steelblue", lw=1.5, alpha=0.9, label="Dados")

    # Região selecionada (se houver)
    if selected_range is not None:
        wl_min, wl_max = selected_range
        ymin, ymax = ax.get_ylim()
        ax.axvspan(wl_min, wl_max, alpha=0.15, color="cyan" if not dark else "yellow", zorder=0)
        ax.axvline(wl_min, color="cyan" if not dark else "yellow", linestyle=":", lw=1, alpha=0.7)
        ax.axvline(wl_max, color="cyan" if not dark else "yellow", linestyle=":", lw=1, alpha=0.7)

    # Curva ajustada (se fornecida)
    if fit_curve is not None and fit_data_plot is not None and fit_wl is not None:
        modelo, curva = fit_curve, fit_data_plot
        color_fit = "yellow" if dark else "red"
        ax.plot(fit_wl, curva, color=color_fit, lw=2, linestyle="--", alpha=0.85, label=f"Ajuste {modelo}")
        ax.legend(loc="upper right", fontsize=8, framealpha=0.8)

    # Remover marcadores antigos
    if hasattr(ax, "markers"):
        for m in ax.markers:
            try:
                m.remove()
            except Exception:
                pass
        ax.markers = []
    else:
        ax.markers = []
    if hasattr(ax, "marker"):
        try:
            ax.marker.remove()
        except Exception:
            pass

    if show_peaks:
        peaks = detectar_picos(spec, prominence=prominence, valley=valley)
        for idx in peaks:
            wl_p = wl_nm[idx]
            int_p = spec_plot[idx]
            marker = ax.scatter(
                wl_p,
                int_p,
                color=color_fg,
                marker=11 if valley else 10,
                zorder=5,
                s=60,
            )
            ax.markers.append(marker)

    ax.set_xlim(wl_nm.min(), wl_nm.max())
    ymin, ymax = np.nanmin(spec_plot), np.nanmax(spec_plot)
    margin = (ymax - ymin) * 0.05 if ymax > ymin else 1.0
    if power_db:
        ax.set_ylim(ymin - margin, ymax + margin)
    else:
        ax.set_ylim(max(0, ymin - margin), ymax + margin)


def main():
    root = tk.Tk()
    root.title("Visualizador de Espectros e Picos")
    root.geometry("900x550")
    root.minsize(700, 450)

    # Dados: lista de (caminho, wl_nm, spec)
    spectra_data = []
    current_index = 0
    prominence = 5.0
    show_peaks = False  # Por padrão picos desabilitados
    show_gradient = False  # Por padrão gradiente desabilitado
    show_power_db = False  # Por padrão intensidade linear (não dB)
    fit_curve_enabled = False  # Ajuste de curva desabilitado por padrão
    fit_model = "gaussian"  # Modelo padrão: "gaussian" ou "lorentzian"
    dark_theme = False
    last_clicked_wl = None  # último pico clicado (para copiar)
    last_clicked_int = None
    fit_info = None  # Informações do ajuste: (modelo, params, r², fwhm)
    selected_range = None  # Região selecionada para ajuste: (wl_min, wl_max) ou None
    span_selector = None  # Widget de seleção de região

    # Figura matplotlib embutida
    fig = Figure(figsize=(8, 4), dpi=100)
    ax = fig.add_subplot(111)
    ax.set_xlabel("Comprimento de onda (nm)")
    ax.set_ylabel("Intensidade (u.a.)")
    ax.grid(True)
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

    # Barra de status / título do arquivo
    status_var = tk.StringVar(value="Nenhum arquivo carregado. Use 'Carregar arquivo(s)'.")
    lbl_status = tk.Label(root, textvariable=status_var, fg="gray", font=("Segoe UI", 10))
    lbl_status.pack(anchor=tk.W, padx=8, pady=(0, 4))

    def carregar_arquivos():
        nonlocal spectra_data, current_index
        paths = filedialog.askopenfilenames(
            title="Selecionar arquivo(s) de espectro",
            filetypes=[("Texto / CSV", "*.txt *.csv"), ("Todos", "*.*")],
        )
        if not paths:
            return
        spectra_data = []
        for path in paths:
            wl_m, spec = ler_dados_arquivo(path)
            if not wl_m or not spec:
                messagebox.showwarning(
                    "Aviso",
                    f"Não foi possível ler dados de:\n{os.path.basename(path)}",
                )
                continue
            wl_nm = np.array(wl_m) * 1e9
            spec = np.array(spec)
            spectra_data.append((path, wl_nm, spec))
        if not spectra_data:
            messagebox.showwarning("Aviso", "Nenhum espectro válido carregado.")
            return
        current_index = 0
        status_var.set(f"Carregados {len(spectra_data)} arquivo(s). Navegue com < e > ou setas.")
        atualizar_grafico()

    def atualizar_grafico():
        nonlocal last_clicked_wl, last_clicked_int, fit_info, span_selector, selected_range
        if not spectra_data:
            ax.clear()
            ax.set_xlabel("Comprimento de onda (nm)")
            ax.set_ylabel("Intensidade (u.a.)")
            ax.text(0.5, 0.5, "Carregue arquivo(s) de espectro.", ha="center", va="center", transform=ax.transAxes)
            canvas.draw_idle()
            return
        idx = max(0, min(current_index, len(spectra_data) - 1))
        path, wl_nm, spec = spectra_data[idx]
        nome = os.path.basename(path)
        
        # Filtra dados se houver região selecionada
        wl_fit = wl_nm
        spec_fit = spec
        if selected_range is not None and fit_curve_enabled:
            wl_min, wl_max = selected_range
            mask = (wl_nm >= wl_min) & (wl_nm <= wl_max)
            if np.sum(mask) > 10:  # Mínimo de pontos para ajuste
                wl_fit = wl_nm[mask]
                spec_fit = spec[mask]
        
        # Ajuste de curva (se habilitado)
        fit_curve_name = None
        fit_curve_data = None
        fit_curve_wl = None
        if fit_curve_enabled:
            params, curva, r2, fwhm = ajustar_curva(wl_fit, spec_fit, modelo=fit_model)
            if params is not None:
                fit_curve_name = "Gaussiana" if fit_model == "gaussian" else "Lorentziana"
                # Curva ajustada sobre a região filtrada
                fit_curve_wl = wl_fit
                fit_curve_data = curva
                fit_info = (fit_model, params, r2, fwhm)
                # Atualiza status com info do ajuste
                amp, center, width = params
                range_info = ""
                if selected_range is not None:
                    range_info = f" [Região: {selected_range[0]:.1f}–{selected_range[1]:.1f}nm]"
                status_var.set(
                    f"Arquivo {idx + 1}/{len(spectra_data)}: {nome}  |  "
                    f"Ajuste {fit_curve_name}: λ={center:.2f}nm, A={amp:.1f}, FWHM={fwhm:.2f}nm, R²={r2:.4f}{range_info}"
                )
            else:
                fit_info = None
                status_var.set(f"Arquivo {idx + 1} / {len(spectra_data)}: {nome}  |  [ERRO] Falha no ajuste de curva")
        else:
            fit_info = None
            status_var.set(f"Arquivo {idx + 1} / {len(spectra_data)}: {nome}")
        
        plotar_espectro_com_picos(
            ax, wl_nm, spec,
            prominence=prominence,
            dark=dark_theme,
            show_peaks=show_peaks,
            show_gradient=show_gradient,
            fit_curve=fit_curve_name,
            fit_data=fit_curve_data,
            fit_wl=fit_curve_wl,
            selected_range=selected_range,
            power_db=show_power_db,
        )
        
        # Ativa/desativa SpanSelector conforme fit_curve_enabled
        if fit_curve_enabled and span_selector is None:
            def onselect(xmin, xmax):
                nonlocal selected_range
                selected_range = (float(xmin), float(xmax))
                atualizar_grafico()
            
            span_selector = SpanSelector(
                ax,
                onselect,
                "horizontal",
                useblit=True,
                props=dict(alpha=0.3, facecolor="cyan"),
                interactive=True,
                drag_from_anywhere=True,
            )
        elif not fit_curve_enabled and span_selector is not None:
            span_selector.set_active(False)
            span_selector = None
            selected_range = None
        
        # Se exibir picos e houver exatamente um, já mostrar suas informações
        if show_peaks and not fit_curve_enabled:
            peaks = detectar_picos(spec, prominence=prominence)
            if len(peaks) == 1:
                wl_p = float(wl_nm[peaks[0]])
                int_p = float(spec[peaks[0]])
                last_clicked_wl = wl_p
                last_clicked_int = int_p
                status_var.set(
                    f"Pico clicado: λ = {wl_p:.2f} nm  |  Intensidade = {int_p:.2f} u.a.  (use os botões para copiar)"
                )
        elif not fit_curve_enabled:
            last_clicked_wl = None
            last_clicked_int = None
        canvas.draw_idle()

    def anterior():
        nonlocal current_index
        if not spectra_data:
            return
        current_index = (current_index - 1) % len(spectra_data)
        atualizar_grafico()

    def proximo():
        nonlocal current_index
        if not spectra_data:
            return
        current_index = (current_index + 1) % len(spectra_data)
        atualizar_grafico()

    def on_key(event):
        key = event.keysym
        if key in ("Left", "Prior", "comma", "less", "minus"):
            anterior()
            return "break"
        if key in ("Right", "Next", "period", "greater", "plus"):
            proximo()
            return "break"

    # Navegação por teclado: funciona em qualquer lugar da janela
    root.bind_all("<Left>", on_key)
    root.bind_all("<Right>", on_key)
    root.bind_all("<Prior>", on_key)    # Page Up
    root.bind_all("<Next>", on_key)     # Page Down
    root.bind_all("<comma>", on_key)    # , (anterior)
    root.bind_all("<period>", on_key)   # . (próximo)
    root.bind_all("<less>", on_key)     # < (anterior)
    root.bind_all("<greater>", on_key)  # > (próximo)
    root.bind_all("<minus>", on_key)    # - (anterior)
    root.bind_all("<plus>", on_key)     # + (próximo)

    def on_click_grafico(event):
        nonlocal last_clicked_wl, last_clicked_int
        if event.inaxes != ax or not spectra_data or not show_peaks:
            return
        idx = max(0, min(current_index, len(spectra_data) - 1))
        _, wl_nm, spec = spectra_data[idx]
        peaks = detectar_picos(spec, prominence=prominence)
        if len(peaks) == 0:
            return
        x_click, y_click = event.xdata, event.ydata
        if x_click is None or y_click is None:
            return
        wl_picos = wl_nm[peaks]
        int_picos = spec[peaks]
        distancias = (wl_picos - x_click) ** 2 + (int_picos - y_click) ** 2
        i_min = int(np.argmin(distancias))
        wl_p = float(wl_picos[i_min])
        int_p = float(int_picos[i_min])
        last_clicked_wl = wl_p
        last_clicked_int = int_p
        status_var.set(
            f"Pico clicado: λ = {wl_p:.2f} nm  |  Intensidade = {int_p:.2f} u.a.  (use os botões para copiar)"
        )

    def copiar_lambda():
        wl_valor = last_clicked_wl
        origem = "pico"
        
        # Se não há pico clicado, tenta usar o centro da curva ajustada
        if wl_valor is None and fit_info is not None:
            modelo, params, r2, fwhm = fit_info
            amp, center, width = params
            wl_valor = center
            origem = "curva ajustada"
        
        if wl_valor is None:
            status_var.set("Clique em um pico ou ajuste uma curva antes de copiar λ.")
            return
        
        texto = f"{wl_valor:.6g}"
        root.clipboard_clear()
        root.clipboard_append(texto)
        root.update()
        status_var.set(f"λ = {texto} nm ({origem}) copiado para a área de transferência.")

    def copiar_intensidade():
        int_valor = last_clicked_int
        origem = "pico"
        
        # Se não há pico clicado, tenta usar a amplitude da curva ajustada
        if int_valor is None and fit_info is not None:
            modelo, params, r2, fwhm = fit_info
            amp, center, width = params
            int_valor = amp
            origem = "curva ajustada"
        
        if int_valor is None:
            status_var.set("Clique em um pico ou ajuste uma curva antes de copiar intensidade.")
            return
        
        texto = f"{int_valor:.6g}"
        root.clipboard_clear()
        root.clipboard_append(texto)
        root.update()
        status_var.set(f"Intensidade = {texto} ({origem}) copiada para a área de transferência.")

    canvas.mpl_connect("button_press_event", on_click_grafico)

    # Botões e controles
    fr_btn = ttk.Frame(root, padding=4)
    fr_btn.pack(fill=tk.X, padx=8, pady=4)
    ttk.Button(fr_btn, text="Carregar arquivo(s)...", command=carregar_arquivos).pack(side=tk.LEFT, padx=(0, 8))
    ttk.Button(fr_btn, text="< Anterior", command=anterior).pack(side=tk.LEFT, padx=2)
    ttk.Button(fr_btn, text="Próximo >", command=proximo).pack(side=tk.LEFT, padx=2)
    #tk.Label(fr_btn, text="Navegação: ← →  |  < >  |  − +  |  Page Up/Down", fg="gray").pack(side=tk.LEFT, padx=12)

    # Exibir picos (por padrão desligado)
    def toggle_peaks():
        nonlocal show_peaks
        show_peaks = var_show_peaks.get()
        atualizar_grafico()

    var_show_peaks = tk.BooleanVar(value=False)
    chk_peaks = ttk.Checkbutton(
        fr_btn,
        text="Exibir picos",
        variable=var_show_peaks,
        command=toggle_peaks,
    )
    chk_peaks.pack(side=tk.LEFT, padx=(16, 4))

    # Gradiente de cores (por padrão desligado)
    def toggle_gradient():
        nonlocal show_gradient
        show_gradient = var_show_gradient.get()
        atualizar_grafico()

    var_show_gradient = tk.BooleanVar(value=False)
    chk_gradient = ttk.Checkbutton(
        fr_btn,
        text="Gradiente de cores",
        variable=var_show_gradient,
        command=toggle_gradient,
    )
    chk_gradient.pack(side=tk.LEFT, padx=(8, 4))

    # Potência em dB (por padrão False)
    def toggle_power_db():
        nonlocal show_power_db
        show_power_db = var_show_power_db.get()
        atualizar_grafico()

    var_show_power_db = tk.BooleanVar(value=False)
    chk_power_db = ttk.Checkbutton(
        fr_btn,
        text="Potência (dB)",
        variable=var_show_power_db,
        command=toggle_power_db,
    )
    chk_power_db.pack(side=tk.LEFT, padx=(8, 4))

    # Ajuste de curva gaussiana/lorentziana
    def toggle_fit():
        nonlocal fit_curve_enabled
        fit_curve_enabled = var_fit_curve.get()
        atualizar_grafico()

    var_fit_curve = tk.BooleanVar(value=False)
    chk_fit = ttk.Checkbutton(
        fr_btn,
        text="Ajustar curva",
        variable=var_fit_curve,
        command=toggle_fit,
    )
    chk_fit.pack(side=tk.LEFT, padx=(8, 4))

    # Modelo de ajuste (Gaussiana ou Lorentziana)
    def change_fit_model(event=None):
        nonlocal fit_model
        fit_model = fit_model_var.get()
        
        if fit_curve_enabled:
            atualizar_grafico()

    fit_model_var = tk.StringVar(value="gaussian")
    combo_model = ttk.Combobox(
        fr_btn,
        textvariable=fit_model_var,
        values=["gaussian", "lorentzian"],
        state="readonly",
        width=10,
    )
    combo_model.pack(side=tk.LEFT, padx=2)
    combo_model.bind("<<ComboboxSelected>>", change_fit_model)

    # Botão para limpar seleção de região
    def limpar_selecao():
        nonlocal selected_range
        selected_range = None
        atualizar_grafico()

    btn_limpar_selecao = ttk.Button(fr_btn, text="Limpar seleção", command=limpar_selecao)
    btn_limpar_selecao.pack(side=tk.LEFT, padx=2)

    # Sensibilidade (prominência): valor maior = menos picos (filtra ruído)
    def on_prominence_change(val=None):
        nonlocal prominence
        if val is None:
            val = prominence_var.get()
        try:
            p = float(val)
            prominence = max(0.5, min(1000.0, p))
            prominence_var.set(prominence)
        except (ValueError, tk.TclError):
            prominence = 5.0
            prominence_var.set(5.0)
        if show_peaks:
            atualizar_grafico()

    tk.Label(fr_btn, text="Prominência (↑ menos picos):", fg="gray").pack(side=tk.LEFT, padx=(8, 2))
    prominence_var = tk.DoubleVar(value=5.0)
    spin_prominence = ttk.Spinbox(
        fr_btn,
        from_=0.5,
        to=200.0,
        increment=0.5,
        textvariable=prominence_var,
        width=6,
        command=on_prominence_change,
    )
    spin_prominence.pack(side=tk.LEFT, padx=2)
    spin_prominence.bind("<Return>", lambda e: on_prominence_change())
    spin_prominence.bind("<FocusOut>", lambda e: on_prominence_change())

    # Copiar λ ou intensidade do último pico clicado
    tk.Label(fr_btn, text="Copiar pico clicado:", fg="gray").pack(side=tk.LEFT, padx=(12, 2))
    ttk.Button(fr_btn, text="Copiar λ", command=copiar_lambda).pack(side=tk.LEFT, padx=2)
    ttk.Button(fr_btn, text="Copiar I", command=copiar_intensidade).pack(side=tk.LEFT, padx=2)

    # Inicial
    atualizar_grafico()
    root.mainloop()


if __name__ == "__main__":
    main()
