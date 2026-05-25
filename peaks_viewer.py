"""
Visualizador de espectros com detecção de picos.

Permite carregar um ou vários arquivos de espectro e visualizar um por vez,
navegando entre eles com < e > (ou setas).

Formatos suportados (auto-detectados):
  - Simples: ``wavelength;intensity`` por linha, sem cabeçalho.
    Wavelength pode estar em metros (~1e-7) ou nm (~100–3000); a detecção é
    automática pela magnitude.
  - ThorLabs FTS (OSA203 e similares): arquivo CSV com bloco ``[SpectrumHeader]``
    contendo metadados ``#Key;Value`` e dados numéricos após ``[Data]``.
    As unidades dos eixos são lidas de ``#XAxisUnit`` (nm_air, nm_vac, m, ...)
    e ``#YAxisUnit`` (dBm, dB, linear, ...).
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
import sys
import json
from pathlib import Path


def _resource_path(rel_path):
    """
    Caminho absoluto para um recurso (ícone, etc.), funcionando tanto rodando
    pelo script quanto a partir de um executável empacotado com PyInstaller
    (``--onefile`` extrai para ``sys._MEIPASS``).
    """
    base = getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base, rel_path)


# ----------------------------------------------------------------------------
#  Persistência de configurações (settings.json)
# ----------------------------------------------------------------------------

def _settings_path():
    """Local do arquivo de configurações (~/.peaksviewer/settings.json)."""
    base = Path.home() / ".peaksviewer"
    try:
        base.mkdir(parents=True, exist_ok=True)
    except Exception:
        # Fallback: ao lado do script/exe se HOME não for gravável
        base = Path(getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__))))
    return base / "settings.json"


def _default_settings():
    """Configurações padrão de fábrica."""
    return {
        "ui_visibility": {
            "navigation": True,
            "peaks": True,
            "gradient": True,
            "power_db": True,
            "fit_curve": True,
            "prominence": True,
            "copy_buttons": True,
        },
        "defaults": {
            "show_peaks": False,
            "show_gradient": False,
            "show_power_db": False,
            "fit_curve_enabled": False,
            "fit_model": "gaussian",
            "auto_enable_db_when_detected": True,
        },
        "appearance": {
            "dark_theme": False,
            "window_geometry": "900x550",
        },
    }


def _merge_settings(loaded, defaults):
    """
    Mescla ``loaded`` em cima de ``defaults`` mantendo apenas chaves conhecidas
    e validando o tipo de cada campo. Garante que um arquivo de settings
    parcialmente desatualizado/quebrado não derrube a aplicação.
    """
    out = {}
    for section, default_section in defaults.items():
        if isinstance(default_section, dict):
            sec = dict(default_section)
            loaded_section = loaded.get(section) if isinstance(loaded, dict) else None
            if isinstance(loaded_section, dict):
                for k, default_v in default_section.items():
                    if k in loaded_section:
                        v = loaded_section[k]
                        if isinstance(default_v, bool) and isinstance(v, bool):
                            sec[k] = v
                        elif isinstance(default_v, (int, float)) and isinstance(v, (int, float)) and not isinstance(v, bool):
                            sec[k] = type(default_v)(v)
                        elif isinstance(default_v, str) and isinstance(v, str):
                            sec[k] = v
            out[section] = sec
        else:
            out[section] = default_section
    return out


def _load_settings():
    """Carrega settings.json (mesclando com os defaults). Nunca lança exceção."""
    path = _settings_path()
    defaults = _default_settings()
    if not path.exists():
        return defaults
    try:
        with open(path, "r", encoding="utf-8") as f:
            return _merge_settings(json.load(f), defaults)
    except Exception:
        return defaults


def _save_settings(settings):
    """Persiste settings.json. Retorna True em sucesso, False caso contrário."""
    try:
        path = _settings_path()
        with open(path, "w", encoding="utf-8") as f:
            json.dump(settings, f, indent=2, ensure_ascii=False)
        return True
    except Exception:
        return False


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
    Lê arquivo de espectro com auto-detecção de formato.

    Suporta o formato simples ``wavelength;intensity`` e o formato ThorLabs FTS
    (cabeçalho com ``#Key;Value`` e seção ``[Data]``). A unidade do eixo X é
    convertida para nm; a unidade do eixo Y é preservada, mas sinalizada via
    ``y_is_db`` quando a fonte estiver em escala logarítmica (dBm/dB/dBW).

    Retorna ``(wl_nm, intensidade, y_is_db)`` ou ``([], [], False)`` em caso de erro.
    """
    wl_raw = []
    int_raw = []
    metadata = {}
    in_data_section = False
    has_data_marker = False

    try:
        with open(caminho_arquivo, "r", encoding="utf-8", errors="replace") as arquivo:
            for linha in arquivo:
                linha = linha.strip()
                if not linha:
                    continue

                # Marcador de seção (ex.: [SpectrumHeader], [Data])
                if linha.startswith("[") and linha.endswith("]"):
                    in_data_section = linha.lower() == "[data]"
                    has_data_marker = has_data_marker or in_data_section
                    continue

                # Linhas de cabeçalho: "#Key;Value"
                if linha.startswith("#"):
                    partes = linha[1:].split(";", 1)
                    if len(partes) == 2:
                        metadata[partes[0].strip()] = partes[1].strip()
                    continue

                # Se o arquivo tem marcador [Data], só lê após ele
                if has_data_marker and not in_data_section:
                    continue

                dados = linha.split(";")
                if len(dados) >= 2:
                    try:
                        wl_raw.append(float(dados[0]))
                        int_raw.append(float(dados[1]))
                    except ValueError:
                        continue
    except FileNotFoundError:
        return [], [], False
    except Exception:
        return [], [], False

    if not wl_raw:
        return [], [], False

    wl_arr = np.array(wl_raw, dtype=float)
    int_arr = np.array(int_raw, dtype=float)

    # Eixo X: garantir nm a partir do cabeçalho ou da magnitude
    x_unit = metadata.get("XAxisUnit", "").lower()
    if x_unit.startswith("nm"):
        wl_nm = wl_arr
    elif x_unit == "m":
        wl_nm = wl_arr * 1e9
    else:
        # Sem dica do cabeçalho: nm fica em ~100–3000; metros em ~1e-7
        if float(np.median(np.abs(wl_arr))) < 1e-3:
            wl_nm = wl_arr * 1e9
        else:
            wl_nm = wl_arr

    # Eixo Y: detecta escala logarítmica (dBm/dB/dBW) pelo cabeçalho
    y_unit = metadata.get("YAxisUnit", "").lower()
    y_is_db = y_unit in ("db", "dbm", "dbw")

    # Heurística para arquivos sem cabeçalho: intensidade linear física é
    # não-negativa, então valores predominantemente negativos indicam dB
    # (ex.: .txt cru exportado de OSAs ThorLabs com dBm).
    if not y_is_db and int_arr.size > 0:
        if np.median(int_arr) < 0 or float(np.mean(int_arr < 0)) > 0.7:
            y_is_db = True

    return wl_nm.tolist(), int_arr.tolist(), y_is_db


def detectar_picos(intensidade, prominence=5, valley=False):
    """Retorna índices dos picos (ou vales) e (wl_nm, intensity) nesses pontos."""
    arr = np.asarray(intensidade)
    if valley:
        peaks, _ = find_peaks(-arr, prominence=prominence)
    else:
        peaks, _ = find_peaks(arr, prominence=prominence)
    return peaks


def _efetivar_spec(spec, y_is_db, power_db):
    """
    Converte o espectro armazenado para a escala efetivamente em uso.

    Single source of truth: o resultado é usado para o plot, detecção de picos,
    ajuste de curva e valores copiados — assim o toggle "Potência (dB)" se
    comporta como um verdadeiro alternador dB ↔ linear sem alterar o arquivo.

    Combinações:
      - y_is_db=True,  power_db=True  → dB rel. ao pico (spec − max(spec))
      - y_is_db=True,  power_db=False → linear, convertido de dB (10^(spec/10))
      - y_is_db=False, power_db=True  → dB rel. ao pico (10·log10(spec/max))
      - y_is_db=False, power_db=False → linear cru (spec)
    """
    arr = np.asarray(spec, dtype=float)
    if arr.size == 0:
        return arr
    if y_is_db:
        if power_db:
            return arr - float(np.max(arr))
        return np.power(10.0, arr / 10.0)
    if power_db:
        ref = float(np.max(arr))
        if ref <= 0:
            ref = 1.0
        return 10.0 * np.log10(np.maximum(arr / ref, 1e-12))
    return arr


def _ylabel_para_escala(y_is_db, power_db):
    """Rótulo do eixo Y para a escala efetiva atual."""
    if power_db:
        return "Potência (dB, rel. ao pico)"
    if y_is_db:
        return "Intensidade linear (convertida de dB)"
    return "Intensidade (u.a.)"


def plotar_espectro_com_picos(ax, wl_nm, spec, prominence=5, valley=False, dark=False, show_peaks=False, show_gradient=False, fit_curve=None, fit_data_plot=None, fit_wl=None, selected_range=None, power_db=False, y_is_db=False):
    """
    Plota o espectro convertido para a escala atualmente em uso (linear ou dB),
    desenhada por ``_efetivar_spec``.

    O argumento ``spec`` é o vetor armazenado tal como foi lido do arquivo;
    a conversão para a escala visível acontece aqui dentro. Já ``fit_data_plot``
    deve vir **já na escala visível** (porque o ajuste de curva é executado
    sobre ``_efetivar_spec(spec, ...)`` em ``atualizar_grafico``).

    Se ``show_peaks=True``, detecta picos no espectro efetivo (mesma escala do
    gráfico, então a prominência tem o mesmo significado do que se vê).
    """
    ax.clear()
    color_fg = "white" if dark else "black"
    color_bg = "black" if dark else "white"
    ax.set_facecolor(color_bg)
    ax.set_xlabel("Comprimento de onda (nm)", color=color_fg)
    ax.tick_params(colors=color_fg)
    ax.grid(True)

    spec_plot = _efetivar_spec(spec, y_is_db, power_db)
    ax.set_ylabel(_ylabel_para_escala(y_is_db, power_db), color=color_fg)
    in_log_scale = bool(power_db)

    if show_gradient:
        gradient_colors = precompute_gradient(wl_nm, dark=dark)
        floor = np.min(spec_plot) - 10 if in_log_scale else 0
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
        peaks = detectar_picos(spec_plot, prominence=prominence, valley=valley)
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
    ymin, ymax = float(np.nanmin(spec_plot)), float(np.nanmax(spec_plot))
    margin = (ymax - ymin) * 0.05 if ymax > ymin else max(abs(ymax) * 0.05, 1e-12)
    if in_log_scale or ymin < 0:
        ax.set_ylim(ymin - margin, ymax + margin)
    else:
        ax.set_ylim(max(0, ymin - margin), ymax + margin)


def _abrir_modal_configuracoes(root, settings, on_apply):
    """
    Modal de Configurações. Edita ``settings`` in-place quando o usuário salva,
    persiste em ``settings.json`` e dispara ``on_apply()`` para o chamador
    aplicar mudanças vivas (visibilidade da toolbar, tema, geometria).

    Estados iniciais (aba "Comportamento") são apenas persistidos — passam a
    valer no próximo arranque, não interferem na sessão atual.
    """
    modal = tk.Toplevel(root)
    modal.title("Configurações")
    modal.transient(root)
    modal.grab_set()
    modal.geometry("560x460")
    modal.minsize(480, 380)

    nb = ttk.Notebook(modal)
    nb.pack(fill=tk.BOTH, expand=True, padx=10, pady=(10, 4))

    # ----- Tab 1: Interface (visibilidade) -----
    tab_iface = ttk.Frame(nb)
    nb.add(tab_iface, text="Interface")

    ttk.Label(
        tab_iface,
        text="Mostrar/ocultar grupos de controles na barra superior:",
        foreground="gray",
    ).pack(anchor=tk.W, padx=12, pady=(12, 4))

    visibility_vars = {}
    for key, label in [
        ("navigation", "Botões de navegação (Anterior / Próximo)"),
        ("peaks", "Checkbox \"Exibir picos\""),
        ("gradient", "Checkbox \"Gradiente de cores\""),
        ("power_db", "Checkbox \"Potência (dB)\""),
        ("fit_curve", "Ajuste de curva (checkbox + modelo + Limpar seleção)"),
        ("prominence", "Spinbox de Prominência"),
        ("copy_buttons", "Botões Copiar λ / Copiar I"),
    ]:
        v = tk.BooleanVar(value=settings["ui_visibility"].get(key, True))
        visibility_vars[key] = v
        ttk.Checkbutton(tab_iface, text=label, variable=v).pack(anchor=tk.W, padx=24, pady=2)

    # ----- Tab 2: Comportamento (defaults) -----
    tab_def = ttk.Frame(nb)
    nb.add(tab_def, text="Comportamento")

    ttk.Label(
        tab_def,
        text="Estados iniciais (aplicados ao abrir o programa):",
        foreground="gray",
    ).pack(anchor=tk.W, padx=12, pady=(12, 4))

    default_vars = {}
    for key, label in [
        ("show_peaks", "Iniciar com \"Exibir picos\" ligado"),
        ("show_gradient", "Iniciar com \"Gradiente de cores\" ligado"),
        ("show_power_db", "Iniciar com \"Potência (dB)\" ligado"),
        ("fit_curve_enabled", "Iniciar com \"Ajustar curva\" ligado"),
    ]:
        v = tk.BooleanVar(value=settings["defaults"].get(key, False))
        default_vars[key] = v
        ttk.Checkbutton(tab_def, text=label, variable=v).pack(anchor=tk.W, padx=24, pady=2)

    fr_model = ttk.Frame(tab_def)
    fr_model.pack(anchor=tk.W, padx=24, pady=(8, 2), fill=tk.X)
    ttk.Label(fr_model, text="Modelo de ajuste padrão:").pack(side=tk.LEFT)
    var_fit_model = tk.StringVar(value=settings["defaults"].get("fit_model", "gaussian"))
    ttk.Combobox(fr_model, textvariable=var_fit_model,
                 values=["gaussian", "lorentzian"], state="readonly", width=14).pack(side=tk.LEFT, padx=(8, 0))

    ttk.Separator(tab_def, orient="horizontal").pack(fill=tk.X, padx=12, pady=10)

    ttk.Label(tab_def, text="Comportamentos automáticos:", foreground="gray").pack(anchor=tk.W, padx=12, pady=(0, 4))
    var_auto_db = tk.BooleanVar(value=settings["defaults"].get("auto_enable_db_when_detected", True))
    ttk.Checkbutton(
        tab_def,
        text="Auto-ativar \"Potência (dB)\" ao carregar dados em escala dB",
        variable=var_auto_db,
    ).pack(anchor=tk.W, padx=24, pady=2)

    # ----- Tab 3: Aparência -----
    tab_apa = ttk.Frame(nb)
    nb.add(tab_apa, text="Aparência")

    var_dark = tk.BooleanVar(value=settings["appearance"].get("dark_theme", False))
    ttk.Checkbutton(tab_apa, text="Tema escuro (gráfico)", variable=var_dark).pack(anchor=tk.W, padx=24, pady=(16, 4))

    fr_geom = ttk.Frame(tab_apa)
    fr_geom.pack(anchor=tk.W, padx=24, pady=(8, 4), fill=tk.X)
    ttk.Label(fr_geom, text="Tamanho inicial da janela (LarguraxAltura):").pack(side=tk.LEFT)
    var_geom = tk.StringVar(value=settings["appearance"].get("window_geometry", "900x550"))
    ttk.Entry(fr_geom, textvariable=var_geom, width=14).pack(side=tk.LEFT, padx=(8, 0))

    ttk.Label(
        tab_apa,
        text=f"Arquivo de configuração: {_settings_path()}",
        foreground="gray",
        wraplength=520,
        justify=tk.LEFT,
    ).pack(anchor=tk.W, padx=12, pady=(20, 4))

    # ----- Botões inferiores -----
    fr_btn = ttk.Frame(modal)
    fr_btn.pack(fill=tk.X, padx=10, pady=(4, 10))

    def restaurar_padroes():
        defaults = _default_settings()
        for k, v in visibility_vars.items():
            v.set(defaults["ui_visibility"][k])
        for k, v in default_vars.items():
            v.set(defaults["defaults"][k])
        var_fit_model.set(defaults["defaults"]["fit_model"])
        var_auto_db.set(defaults["defaults"]["auto_enable_db_when_detected"])
        var_dark.set(defaults["appearance"]["dark_theme"])
        var_geom.set(defaults["appearance"]["window_geometry"])

    def salvar_e_aplicar():
        # Atualiza ``settings`` in-place
        for k, v in visibility_vars.items():
            settings["ui_visibility"][k] = bool(v.get())
        for k, v in default_vars.items():
            settings["defaults"][k] = bool(v.get())
        settings["defaults"]["fit_model"] = var_fit_model.get()
        settings["defaults"]["auto_enable_db_when_detected"] = bool(var_auto_db.get())
        settings["appearance"]["dark_theme"] = bool(var_dark.get())
        # Validação simples para a geometria (ex.: "900x550")
        geom_txt = var_geom.get().strip()
        if "x" in geom_txt:
            try:
                w, h = geom_txt.lower().split("x", 1)
                int(w); int(h)
                settings["appearance"]["window_geometry"] = f"{int(w)}x{int(h)}"
            except ValueError:
                pass

        # Aplica geometria imediatamente (apenas largura/altura — nunca move a janela)
        try:
            root.geometry(settings["appearance"]["window_geometry"])
        except Exception:
            pass

        ok = _save_settings(settings)
        if not ok:
            messagebox.showwarning(
                "Configurações",
                "Não foi possível gravar settings.json — alterações aplicadas apenas nesta sessão.",
            )

        on_apply()
        modal.destroy()

    ttk.Button(fr_btn, text="Restaurar padrões", command=restaurar_padroes).pack(side=tk.LEFT)
    ttk.Button(fr_btn, text="Cancelar", command=modal.destroy).pack(side=tk.RIGHT)
    ttk.Button(fr_btn, text="Salvar e aplicar", command=salvar_e_aplicar).pack(side=tk.RIGHT, padx=(0, 6))

    modal.bind("<Escape>", lambda e: modal.destroy())
    modal.wait_window()


def main():
    settings = _load_settings()

    root = tk.Tk()
    root.title("Visualizador de Espectros e Picos")
    root.geometry(settings["appearance"].get("window_geometry", "900x550"))
    root.minsize(700, 450)

    # Ícone da janela: tenta .ico (Windows nativo), depois .png como fallback.
    # Silencioso se ambos falharem — o app continua funcionando sem ícone custom.
    try:
        ico = _resource_path(os.path.join("assets", "logo.ico"))
        if os.path.exists(ico):
            root.iconbitmap(default=ico)
        else:
            png = _resource_path(os.path.join("assets", "logo.png"))
            if os.path.exists(png):
                root.iconphoto(True, tk.PhotoImage(file=png))
    except Exception:
        pass

    # Dados: lista de (caminho, wl_nm, spec, y_is_db)
    spectra_data = []
    current_index = 0
    prominence = 5.0
    last_spec_range = 0.0  # range do spec_eff na última recalculação (escala atual)
    # Estados iniciais vêm de settings["defaults"]; o usuário ajusta via toolbar
    # e configura novos defaults pelo modal "Configurações".
    show_peaks = bool(settings["defaults"].get("show_peaks", False))
    show_gradient = bool(settings["defaults"].get("show_gradient", False))
    show_power_db = bool(settings["defaults"].get("show_power_db", False))
    fit_curve_enabled = bool(settings["defaults"].get("fit_curve_enabled", False))
    fit_model = str(settings["defaults"].get("fit_model", "gaussian"))
    dark_theme = bool(settings["appearance"].get("dark_theme", False))
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
        nonlocal spectra_data, current_index, show_power_db
        paths = filedialog.askopenfilenames(
            title="Selecionar arquivo(s) de espectro",
            filetypes=[
                ("Texto / CSV", "*.txt *.csv"),
                ("ThorLabs FTS", "*.csv *.spf2 *.txt"),
                ("Todos", "*.*"),
            ],
        )
        if not paths:
            return
        spectra_data = []
        for path in paths:
            wl_list, spec, y_is_db = ler_dados_arquivo(path)
            if not wl_list or len(spec) == 0:
                messagebox.showwarning(
                    "Aviso",
                    f"Não foi possível ler dados de:\n{os.path.basename(path)}",
                )
                continue
            wl_nm = np.array(wl_list)
            spec = np.array(spec)
            spectra_data.append((path, wl_nm, spec, y_is_db))
        if not spectra_data:
            messagebox.showwarning("Aviso", "Nenhum espectro válido carregado.")
            return
        current_index = 0

        # Auto-ativa "Potência (dB)" quando algum arquivo já vem em escala log
        # (ex.: ThorLabs FTS com #YAxisUnit;dBm). Pode ser desabilitado em
        # Configurações → Comportamento.
        if (
            settings["defaults"].get("auto_enable_db_when_detected", True)
            and any(s[3] for s in spectra_data)
            and not show_power_db
        ):
            show_power_db = True
            var_show_power_db.set(True)

        # Reset da prominência default e da faixa do spinbox para a nova escala
        _recalc_prominence(force_default=True)

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
        path, wl_nm, spec, y_is_db = spectra_data[idx]
        nome = os.path.basename(path)

        # Escala efetiva (linear ou dB) usada por TODOS os cálculos a partir daqui:
        # plot, picos, ajuste e valores copiados. Quando o usuário alterna o toggle,
        # tudo recalcula automaticamente sobre a nova escala.
        spec_eff = _efetivar_spec(spec, y_is_db, show_power_db)

        # Filtra dados se houver região selecionada (ajuste é feito em spec_eff)
        wl_fit = wl_nm
        spec_fit = spec_eff
        if selected_range is not None and fit_curve_enabled:
            wl_min, wl_max = selected_range
            mask = (wl_nm >= wl_min) & (wl_nm <= wl_max)
            if np.sum(mask) > 10:  # Mínimo de pontos para ajuste
                wl_fit = wl_nm[mask]
                spec_fit = spec_eff[mask]

        # Unidade exibida ao usuário (para textos de status / botões de copiar)
        unidade = "dB" if show_power_db else ("u.a." if not y_is_db else "lin")

        # Ajuste de curva (se habilitado) — opera sobre spec_eff
        fit_curve_name = None
        fit_curve_data = None
        fit_curve_wl = None
        if fit_curve_enabled:
            params, curva, r2, fwhm = ajustar_curva(wl_fit, spec_fit, modelo=fit_model)
            if params is not None:
                fit_curve_name = "Gaussiana" if fit_model == "gaussian" else "Lorentziana"
                fit_curve_wl = wl_fit
                fit_curve_data = curva  # já está na escala efetiva
                fit_info = (fit_model, params, r2, fwhm)
                amp, center, width = params
                range_info = ""
                if selected_range is not None:
                    range_info = f" [Região: {selected_range[0]:.1f}–{selected_range[1]:.1f}nm]"
                status_var.set(
                    f"Arquivo {idx + 1}/{len(spectra_data)}: {nome}  |  "
                    f"Ajuste {fit_curve_name}: λ={center:.2f}nm, A={amp:.4g} {unidade}, "
                    f"FWHM={fwhm:.2f}nm, R²={r2:.4f}{range_info}"
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
            fit_data_plot=fit_curve_data,
            fit_wl=fit_curve_wl,
            selected_range=selected_range,
            power_db=show_power_db,
            y_is_db=y_is_db,
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
            peaks = detectar_picos(spec_eff, prominence=prominence)
            if len(peaks) == 1:
                wl_p = float(wl_nm[peaks[0]])
                int_p = float(spec_eff[peaks[0]])
                last_clicked_wl = wl_p
                last_clicked_int = int_p
                status_var.set(
                    f"Pico clicado: λ = {wl_p:.2f} nm  |  Intensidade = {int_p:.4g} {unidade}  (use os botões para copiar)"
                )
        elif not fit_curve_enabled:
            last_clicked_wl = None
            last_clicked_int = None
        canvas.draw_idle()

    def _recalc_prominence(force_default=False):
        """
        Ajusta prominência e faixa do spinbox para a escala efetiva atual.

        - ``force_default=True``: define prominência como 2% do range do espectro
          efetivo (usado ao carregar arquivos novos).
        - ``force_default=False``: escala a prominência atual proporcionalmente
          ao novo range (usado quando o usuário alterna "Potência (dB)"),
          preservando a sensibilidade relativa.
        """
        nonlocal prominence, last_spec_range
        if not spectra_data:
            return
        idx = max(0, min(current_index, len(spectra_data) - 1))
        _, _, spec_raw, y_is_db_i = spectra_data[idx]
        spec_eff = _efetivar_spec(spec_raw, y_is_db_i, show_power_db)
        rng = float(np.nanmax(spec_eff) - np.nanmin(spec_eff))
        if rng <= 0:
            return

        if force_default or last_spec_range <= 0:
            new_prom = rng * 0.02
        else:
            new_prom = prominence * (rng / last_spec_range)

        # Limites razoáveis para o spinbox
        minp = rng * 1e-5
        maxp = rng * 0.5
        new_prom = max(minp, min(maxp, new_prom))

        prominence = new_prom
        last_spec_range = rng
        prominence_var.set(prominence)
        spin_prominence.configure(from_=minp, to=maxp, increment=rng * 0.005)

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
        _, wl_nm, spec, y_is_db_i = spectra_data[idx]
        spec_eff = _efetivar_spec(spec, y_is_db_i, show_power_db)
        peaks = detectar_picos(spec_eff, prominence=prominence)
        if len(peaks) == 0:
            return
        x_click, y_click = event.xdata, event.ydata
        if x_click is None or y_click is None:
            return
        wl_picos = wl_nm[peaks]
        int_picos = spec_eff[peaks]
        distancias = (wl_picos - x_click) ** 2 + (int_picos - y_click) ** 2
        i_min = int(np.argmin(distancias))
        wl_p = float(wl_picos[i_min])
        int_p = float(int_picos[i_min])
        last_clicked_wl = wl_p
        last_clicked_int = int_p
        unidade = "dB" if show_power_db else ("u.a." if not y_is_db_i else "lin")
        status_var.set(
            f"Pico clicado: λ = {wl_p:.2f} nm  |  Intensidade = {int_p:.4g} {unidade}  (use os botões para copiar)"
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

    # Toolbar é construída como uma lista ordenada de widgets para permitir
    # mostrar/ocultar grupos a partir do modal Configurações sem perder a ordem.
    toolbar_layout = []  # lista de tuplas: (group_key|None, widget, pack_kwargs)

    btn_carregar = ttk.Button(fr_btn, text="Carregar arquivo(s)...", command=carregar_arquivos)
    toolbar_layout.append((None, btn_carregar, dict(side=tk.LEFT, padx=(0, 8))))

    btn_anterior = ttk.Button(fr_btn, text="< Anterior", command=anterior)
    btn_proximo = ttk.Button(fr_btn, text="Próximo >", command=proximo)
    toolbar_layout.append(("navigation", btn_anterior, dict(side=tk.LEFT, padx=2)))
    toolbar_layout.append(("navigation", btn_proximo, dict(side=tk.LEFT, padx=2)))

    # Exibir picos
    def toggle_peaks():
        nonlocal show_peaks
        show_peaks = var_show_peaks.get()
        atualizar_grafico()

    var_show_peaks = tk.BooleanVar(value=show_peaks)
    chk_peaks = ttk.Checkbutton(fr_btn, text="Exibir picos", variable=var_show_peaks, command=toggle_peaks)
    toolbar_layout.append(("peaks", chk_peaks, dict(side=tk.LEFT, padx=(16, 4))))

    # Gradiente de cores
    def toggle_gradient():
        nonlocal show_gradient
        show_gradient = var_show_gradient.get()
        atualizar_grafico()

    var_show_gradient = tk.BooleanVar(value=show_gradient)
    chk_gradient = ttk.Checkbutton(fr_btn, text="Gradiente de cores", variable=var_show_gradient, command=toggle_gradient)
    toolbar_layout.append(("gradient", chk_gradient, dict(side=tk.LEFT, padx=(8, 4))))

    # Potência em dB
    def toggle_power_db():
        nonlocal show_power_db
        show_power_db = var_show_power_db.get()
        _recalc_prominence(force_default=False)
        atualizar_grafico()

    var_show_power_db = tk.BooleanVar(value=show_power_db)
    chk_power_db = ttk.Checkbutton(fr_btn, text="Potência (dB)", variable=var_show_power_db, command=toggle_power_db)
    toolbar_layout.append(("power_db", chk_power_db, dict(side=tk.LEFT, padx=(8, 4))))

    # Ajuste de curva gaussiana/lorentziana + combobox + limpar seleção
    def toggle_fit():
        nonlocal fit_curve_enabled
        fit_curve_enabled = var_fit_curve.get()
        atualizar_grafico()

    var_fit_curve = tk.BooleanVar(value=fit_curve_enabled)
    chk_fit = ttk.Checkbutton(fr_btn, text="Ajustar curva", variable=var_fit_curve, command=toggle_fit)
    toolbar_layout.append(("fit_curve", chk_fit, dict(side=tk.LEFT, padx=(8, 4))))

    def change_fit_model(event=None):
        nonlocal fit_model
        fit_model = fit_model_var.get()
        if fit_curve_enabled:
            atualizar_grafico()

    fit_model_var = tk.StringVar(value=fit_model)
    combo_model = ttk.Combobox(fr_btn, textvariable=fit_model_var,
                                values=["gaussian", "lorentzian"], state="readonly", width=10)
    combo_model.bind("<<ComboboxSelected>>", change_fit_model)
    toolbar_layout.append(("fit_curve", combo_model, dict(side=tk.LEFT, padx=2)))

    def limpar_selecao():
        nonlocal selected_range
        selected_range = None
        atualizar_grafico()

    btn_limpar_selecao = ttk.Button(fr_btn, text="Limpar seleção", command=limpar_selecao)
    toolbar_layout.append(("fit_curve", btn_limpar_selecao, dict(side=tk.LEFT, padx=2)))

    # Sensibilidade (prominência)
    def on_prominence_change(val=None):
        nonlocal prominence
        if val is None:
            val = prominence_var.get()
        try:
            p = float(val)
            prominence = max(0.0, p)
            prominence_var.set(prominence)
        except (ValueError, tk.TclError):
            _recalc_prominence(force_default=True)
            return
        if show_peaks:
            atualizar_grafico()

    lbl_prominence = tk.Label(fr_btn, text="Prominência:", fg="gray")
    prominence_var = tk.DoubleVar(value=prominence)
    spin_prominence = ttk.Spinbox(fr_btn, from_=0.5, to=200.0, increment=0.5,
                                   textvariable=prominence_var, width=6,
                                   command=on_prominence_change)
    spin_prominence.bind("<Return>", lambda e: on_prominence_change())
    spin_prominence.bind("<FocusOut>", lambda e: on_prominence_change())
    toolbar_layout.append(("prominence", lbl_prominence, dict(side=tk.LEFT, padx=(8, 2))))
    toolbar_layout.append(("prominence", spin_prominence, dict(side=tk.LEFT, padx=2)))

    # Copiar λ / I
    lbl_copiar = tk.Label(fr_btn, text="Pico:", fg="gray")
    btn_copiar_lambda = ttk.Button(fr_btn, text="Copiar λ", command=copiar_lambda)
    btn_copiar_int = ttk.Button(fr_btn, text="Copiar I", command=copiar_intensidade)
    toolbar_layout.append(("copy_buttons", lbl_copiar, dict(side=tk.LEFT, padx=(12, 2))))
    toolbar_layout.append(("copy_buttons", btn_copiar_lambda, dict(side=tk.LEFT, padx=2)))
    toolbar_layout.append(("copy_buttons", btn_copiar_int, dict(side=tk.LEFT, padx=2)))

    def aplicar_visibilidade_toolbar():
        """Repacka a toolbar respeitando settings['ui_visibility']."""
        for _grp, w, _kw in toolbar_layout:
            try:
                w.pack_forget()
            except Exception:
                pass
        try:
            btn_configuracoes.pack_forget()
        except Exception:
            pass
        # Botão de Configurações sempre visível, ancorado à direita
        btn_configuracoes.pack(side=tk.RIGHT, padx=(8, 0))
        for grp, w, kwargs in toolbar_layout:
            if grp is None or settings["ui_visibility"].get(grp, True):
                w.pack(**kwargs)

    def abrir_configuracoes():
        _abrir_modal_configuracoes(
            root, settings,
            on_apply=lambda: (
                aplicar_settings_runtime(),
                aplicar_visibilidade_toolbar(),
                atualizar_grafico(),
            ),
        )

    def aplicar_settings_runtime():
        """Aplica em tempo real os settings que afetam o estado vivo."""
        nonlocal dark_theme
        dark_theme = bool(settings["appearance"].get("dark_theme", False))

    btn_configuracoes = ttk.Button(fr_btn, text="⚙ Configurações", command=abrir_configuracoes)

    # Aplica visibilidade inicial conforme settings.json
    aplicar_visibilidade_toolbar()

    # Inicial
    atualizar_grafico()
    root.mainloop()


if __name__ == "__main__":
    main()
