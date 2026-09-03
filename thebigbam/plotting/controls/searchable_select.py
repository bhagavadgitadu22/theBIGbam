import json

import param
from panel.custom import JSComponent


def decode_search_request(payload: str) -> tuple[int, str]:
    """Decode the browser's atomic autocomplete request."""
    try:
        request = json.loads(payload)
        nonce = int(request["nonce"])
        query = str(request["query"])
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        raise ValueError("invalid autocomplete request") from error
    return nonce, query


class SearchableSelect(JSComponent):
    """Searchable select dropdown using Tom Select."""

    value = param.String(default="")
    options = param.List(default=[])
    placeholder = param.String(default="")
    server_search = param.Boolean(default=False)
    allow_custom = param.Boolean(default=False)
    search_query = param.String(default="")
    # Query and sequence travel in one property update so the server can never
    # observe the nonce from one request with the text from another.
    search_request = param.String(default="")
    search_result_nonce = param.Integer(default=0)
    search_result_query = param.String(default="")
    scope_nonce = param.Integer(default=0)
    min_search_chars = param.Integer(default=0, bounds=(0, None))
    disabled = param.Boolean(default=False)

    _stylesheets = [
        "https://cdn.jsdelivr.net/npm/tom-select@2.4.1/dist/css/tom-select.css",
        # Override styles to match Bokeh widget look
        """
        .ts-wrapper { width: 100%; }
        .ts-wrapper .ts-control { border: 1px solid #ccc; border-radius: 4px; padding: 4px 8px; min-height: 31px; font-size: 14px; }
        .ts-dropdown { margin-top: 0 !important; }
        .ts-dropdown .ts-dropdown-content { max-height: 350px; }
        """,
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
        if (!document.getElementById('thebigbam-tom-select-css')) {
            const link = document.createElement('link');
            link.id = 'thebigbam-tom-select-css';
            link.rel = 'stylesheet';
            link.href = 'https://cdn.jsdelivr.net/npm/tom-select@2.4.1/dist/css/tom-select.css';
            document.head.appendChild(link);
        }
        if (!document.getElementById('thebigbam-tom-select-overlay-css')) {
            const style = document.createElement('style');
            style.id = 'thebigbam-tom-select-overlay-css';
            style.textContent = `
                body > .ts-dropdown { z-index: 10000; }
                body > .ts-dropdown .ts-dropdown-content { max-height: 350px; }
                .thebigbam-autocomplete-loading {
                    align-items: center; color: #666; display: flex; gap: 8px; padding: 8px;
                }
                .thebigbam-autocomplete-spinner {
                    animation: thebigbam-autocomplete-spin .7s linear infinite;
                    border: 2px solid #ddd; border-radius: 50%; border-top-color: #fc62b8;
                    display: inline-block; height: 12px; width: 12px;
                }
                @keyframes thebigbam-autocomplete-spin { to { transform: rotate(360deg); } }
            `;
            document.head.appendChild(style);
        }
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

        // When server_search is on, typed queries are sent to Python (via
        // model.search_query) instead of being filtered purely client-side —
        // the 'options' handler below resolves the pending load() callback
        // once Python pushes back the matching results.
        const pendingFilterLoads = new Map();
        const queryResults = new Map();
        let requestSequence = model.search_result_nonce;
        let pendingInputQuery = null;
        queryResults.set(model.search_query || '', model.options.slice());
        const tsConfig = {
            create: model.allow_custom,
            persist: false,
            maxOptions: 100,
            placeholder: model.placeholder,
            options: allOptions,
            items: model.value ? [model.value] : [],
            onChange: (val) => { model.value = val; },
            // The coloring-rule list is deliberately scrollable. Portalling
            // the menu to body prevents that ancestor from clipping it.
            dropdownParent: 'body'
        };
        if (model.server_search) {
            tsConfig.loadThrottle = 300;
            tsConfig.shouldLoad = (query) => query.length >= model.min_search_chars;
            tsConfig.load = (query, callback) => {
                const cached = queryResults.get(query);
                if (cached) {
                    pendingInputQuery = null;
                    replaceOptions(cached);
                    callback(cached.map(o => ({value: o, text: o})));
                    return;
                }
                const previous = pendingFilterLoads.get(query);
                if (previous) previous.callback([]);
                model.search_query = query;
                requestSequence += 1;
                pendingFilterLoads.set(query, {callback, nonce: requestSequence});
                model.search_request = JSON.stringify({nonce: requestSequence, query});
            };
            tsConfig.onType = (query) => {
                // Substring operators filter by the literal text, whether or
                // not the user chooses one of the suggested complete values.
                if (model.allow_custom) model.value = query;
                if (query.length >= model.min_search_chars) {
                    const cached = queryResults.get(query);
                    if (cached) {
                        pendingInputQuery = null;
                        replaceOptions(cached);
                    } else {
                        pendingInputQuery = query;
                        replaceOptions([], false);
                        ts.refreshOptions(true);
                    }
                    return;
                }
                pendingInputQuery = null;
                cancelPendingLoads();
                replaceOptions([], false);
            };
            tsConfig.render = {
                no_results: (data, escape) => {
                    const query = data.input || '';
                    if (pendingInputQuery === query) {
                        return '<div class="thebigbam-autocomplete-loading">' +
                            '<span class="thebigbam-autocomplete-spinner"></span>Loading…</div>';
                    }
                    const message = query.length < model.min_search_chars
                        ? `Type at least ${model.min_search_chars} characters`
                        : 'No results found';
                    return `<div class="no-results">${escape(message)}</div>`;
                },
            };
        }

        const ts = new TomSelect(select, tsConfig);
        if (model.disabled) ts.disable();

        function replaceOptions(newOptions, includeValue = true) {
            const newSet = new Set(newOptions);
            // The active value is model state, not part of the bounded search
            // result window. Keep it renderable without growing model.options
            // beyond the server's result limit.
            if (includeValue && model.value) newSet.add(model.value);
            Object.keys(ts.options).forEach((key) => {
                if (!newSet.has(key)) ts.removeOption(key, true);
            });
            newOptions.forEach((opt) => {
                if (!ts.options.hasOwnProperty(opt)) ts.addOption({value: opt, text: opt});
            });
            if (includeValue && model.value && !ts.options.hasOwnProperty(model.value)) {
                ts.addOption({value: model.value, text: model.value});
            }
            allOptions = newOptions.map(o => ({value: o, text: o}));
            ts.refreshOptions(false);
        }

        function cancelPendingLoads() {
            pendingFilterLoads.forEach((request) => request.callback([]));
            pendingFilterLoads.clear();
        }

        model.on('disabled', () => {
            if (model.disabled) ts.disable();
            else ts.enable();
        });

        model.on('allow_custom', () => {
            ts.settings.create = model.allow_custom;
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
            // Search responses are associated with their exact query by the
            // nonce event below. Scope/programmatic replacements have no
            // pending load and can be projected immediately.
            if (!model.server_search || pendingFilterLoads.size === 0) {
                replaceOptions(model.options);
            }
        });

        model.on('search_result_nonce', () => {
            const query = model.search_result_query;
            const currentQuery = ts.control_input.value;
            const request = pendingFilterLoads.get(query);
            if (request && request.nonce === model.search_result_nonce) pendingFilterLoads.delete(query);
            const newOptions = model.options.slice();
            if (!request || request.nonce !== model.search_result_nonce) return;
            queryResults.set(query, newOptions);
            if (query !== currentQuery || query.length < model.min_search_chars) {
                request.callback([]);
                return;
            }
            pendingInputQuery = null;
            replaceOptions(newOptions);
            request.callback(newOptions.map(o => ({value: o, text: o})));
            for (const [oldQuery, oldRequest] of pendingFilterLoads.entries()) {
                if (oldQuery !== currentQuery) {
                    oldRequest.callback([]);
                    pendingFilterLoads.delete(oldQuery);
                }
            }
        });

        model.on('scope_nonce', () => {
            // A MAG/contig/filter/view change invalidates both the visible
            // query and every pending response from the previous scope.
            cancelPendingLoads();
            queryResults.clear();
            if (typeof ts.clearCache === 'function') ts.clearCache();
            ts.setTextboxValue('');
            replaceOptions(model.options);
            queryResults.set('', model.options.slice());
        });

        model.on('value', () => {
            // Echoing each free-text keystroke back through setValue would
            // turn it into a selected item and collapse the text editor.
            if (model.allow_custom && ts.control_input.matches(':focus')) return;
            // Ensure the option exists before setting value
            if (model.value && !ts.options[model.value]) {
                ts.addOption({value: model.value, text: model.value});
            }
            ts.setValue(model.value, true);
        });

        return container;
    }
    """
