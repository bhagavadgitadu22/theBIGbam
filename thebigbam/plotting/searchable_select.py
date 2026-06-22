import param
from panel.custom import JSComponent


class SearchableSelect(JSComponent):
    """Searchable select dropdown using Tom Select."""

    value = param.String(default="")
    options = param.List(default=[])
    placeholder = param.String(default="")

    _stylesheets = [
        "https://cdn.jsdelivr.net/npm/tom-select@2.4.1/dist/css/tom-select.css",
        # Override styles to match Bokeh widget look
        """
        .ts-wrapper { width: 100%; }
        .ts-wrapper .ts-control { border: 1px solid #ccc; border-radius: 4px; padding: 4px 8px; min-height: 31px; font-size: 14px; }
        .ts-dropdown { margin-top: 0 !important; }
        .ts-dropdown .ts-dropdown-content { max-height: 350px; }
        """
    ]

    _esm = """
    import TomSelect from "https://cdn.jsdelivr.net/npm/tom-select@2.4.1/+esm";

    /** Copy text to the system clipboard (best-effort). */
    function copyToClipboard(text) {
        if (navigator.clipboard && navigator.clipboard.writeText) {
            navigator.clipboard.writeText(text).catch(() => {});
        } else {
            // Fallback for older browsers / non-HTTPS contexts
            const ta = document.createElement('textarea');
            ta.value = text;
            ta.style.position = 'fixed';
            ta.style.opacity = '0';
            document.body.appendChild(ta);
            ta.focus();
            ta.select();
            try { document.execCommand('copy'); } catch (_) {}
            document.body.removeChild(ta);
        }
    }

    function showCopiedTooltip(anchorEl) {
        const rect = anchorEl.getBoundingClientRect();
        const tip = document.createElement('div');
        tip.textContent = 'Copied';
        Object.assign(tip.style, {
            position: 'fixed',
            left: (rect.left + rect.width / 2) + 'px',
            top: (rect.top - 28) + 'px',
            transform: 'translateX(-50%)',
            background: '#333',
            color: '#fff',
            padding: '3px 8px',
            borderRadius: '4px',
            fontSize: '12px',
            pointerEvents: 'none',
            zIndex: '9999',
            opacity: '1',
            transition: 'opacity 0.4s ease',
        });
        document.body.appendChild(tip);
        setTimeout(() => { tip.style.opacity = '0'; }, 700);
        setTimeout(() => { if (tip.parentNode) tip.parentNode.removeChild(tip); }, 1100);
    }

    export function render({ model }) {
        const container = document.createElement('div');
        container.style.width = '100%';
        const select = document.createElement('select');
        select.setAttribute('placeholder', model.placeholder);
        container.appendChild(select);

        let allOptions = model.options.map(o => ({value: o, text: o}));

        // Seed the pre-set value so it displays correctly
        const initialOptions = [];
        if (model.value) {
            initialOptions.push({value: model.value, text: model.value});
        }

        const ts = new TomSelect(select, {
            create: false,
            maxOptions: 100,
            placeholder: model.placeholder,
            options: allOptions,
            items: model.value ? [model.value] : [],
            onChange: (val) => { model.value = val; }
        });

        let lastClickTime = 0;
        container.addEventListener('mousedown', (event) => {
            if (event.button !== 0) return;
            const control = container.querySelector('.ts-control');
            if (!control || !control.contains(event.target)) return;
            const now = Date.now();
            if (now - lastClickTime < 400) {
                event.preventDefault();
                event.stopImmediatePropagation();
                const val = model.value;
                if (val) {
                    copyToClipboard(val);
                    if (ts.isOpen) ts.close();
                    showCopiedTooltip(control);
                }
                lastClickTime = 0;
                return;
            }
            lastClickTime = now;
        });

        model.on('options', () => {
            const currentVal = model.value;
            allOptions = model.options.map(o => ({value: o, text: o}));
            // Rebuild options: clear render cache, repopulate, then sync UI
            ts.clearOptions();
            ts.addOptions(allOptions);
            ts.refreshItems();
            ts.refreshOptions(false);  // false = don't open dropdown
            // Restore the current value if still valid
            if (currentVal && model.options.includes(currentVal)) {
                ts.setValue(currentVal, true);
            } else {
                ts.setValue('', true);
            }
        });

        model.on('value', () => {
            // Ensure the option exists before setting value
            if (model.value && !ts.options[model.value]) {
                ts.addOption({value: model.value, text: model.value});
            }
            ts.setValue(model.value, true);
        });

        return container;
    }
    """
