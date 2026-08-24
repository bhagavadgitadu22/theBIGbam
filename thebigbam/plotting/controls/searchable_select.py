import param
from panel.custom import JSComponent


class SearchableSelect(JSComponent):
    """Searchable select dropdown using Tom Select."""

    value = param.String(default="")
    options = param.List(default=[])
    placeholder = param.String(default="")
    server_search = param.Boolean(default=False)
    search_query = param.String(default="")
    search_nonce = param.Integer(default=0)
    search_result_nonce = param.Integer(default=0)
    search_result_query = param.String(default="")
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
        let pendingLoadCallback = null;
        const pendingFilterLoads = new Map();
        const tsConfig = {
            create: false,
            maxOptions: 100,
            placeholder: model.placeholder,
            options: allOptions,
            items: model.value ? [model.value] : [],
            onChange: (val) => { model.value = val; }
        };
        if (model.server_search) {
            tsConfig.loadThrottle = 300;
            tsConfig.shouldLoad = (query) => query.length >= model.min_search_chars;
            tsConfig.load = (query, callback) => {
                if (model.min_search_chars > 0) {
                    const previous = pendingFilterLoads.get(query);
                    if (previous) previous([]);
                    pendingFilterLoads.set(query, callback);
                }
                else pendingLoadCallback = callback;
                model.search_query = query;
                // search_nonce always changes, so Python's watcher fires even
                // when the same text is searched twice in a row (Param skips
                // no-op assignments to search_query itself).
                model.search_nonce = model.search_nonce + 1;
            };
            tsConfig.onType = (query) => {
                if (query.length >= model.min_search_chars) return;
                pendingFilterLoads.forEach((callback) => callback([]));
                pendingFilterLoads.clear();
                ts.clearOptions();
                ts.refreshOptions(false);
            };
            if (model.min_search_chars > 0) {
                tsConfig.render = {
                    no_results: (data, escape) => {
                        const query = data.input || '';
                        const message = query.length < model.min_search_chars
                            ? `Type at least ${model.min_search_chars} characters`
                            : 'No results found';
                        return `<div class="no-results">${escape(message)}</div>`;
                    },
                };
            }
        }

        const ts = new TomSelect(select, tsConfig);
        if (model.disabled) ts.disable();

        model.on('disabled', () => {
            if (model.disabled) ts.disable();
            else ts.enable();
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
            const newOptions = model.options;  // fresh array of strings from Python
            allOptions = newOptions.map(o => ({value: o, text: o}));

            // Surgically diff the option pool only — never touch selection
            // here. model.on('value', ...) below is the sole owner of what's
            // selected; Python already decided whether the value needs to
            // change as part of the same update, so duplicating that
            // decision here (via clear()/clearOptions()/setValue()) is both
            // unnecessary and a source of subtle selection-clobbering bugs.
            const newSet = new Set(newOptions);
            Object.keys(ts.options).forEach((key) => {
                if (!newSet.has(key)) {
                    ts.removeOption(key, true);  // silent: never fires onChange
                }
            });
            newOptions.forEach((opt) => {
                if (!ts.options.hasOwnProperty(opt)) {
                    ts.addOption({value: opt, text: opt});
                }
            });
            ts.refreshOptions(false);  // false = don't open dropdown

            // Resolve any in-flight server-side search: tells Tom Select the
            // load() for the current query has finished, so it renders these
            // results (and clears its loading indicator).
            if (model.min_search_chars <= 0 && pendingLoadCallback) {
                const cb = pendingLoadCallback;
                pendingLoadCallback = null;
                cb(allOptions);
            }
        });

        model.on('search_result_nonce', () => {
            if (model.min_search_chars <= 0) return;
            const query = model.search_result_query;
            const currentQuery = ts.control_input.value;
            const callback = pendingFilterLoads.get(query);
            pendingFilterLoads.delete(query);
            if (!callback) return;
            if (query !== currentQuery || query.length < model.min_search_chars) {
                callback([]);
                return;
            }
            const newOptions = model.options;
            const newSet = new Set(newOptions);
            Object.keys(ts.options).forEach((key) => {
                if (!newSet.has(key)) ts.removeOption(key, true);
            });
            newOptions.forEach((opt) => {
                if (!ts.options.hasOwnProperty(opt)) ts.addOption({value: opt, text: opt});
            });
            callback(newOptions.map(o => ({value: o, text: o})));
            for (const [oldQuery, oldCallback] of pendingFilterLoads.entries()) {
                if (oldQuery !== currentQuery) {
                    oldCallback([]);
                    pendingFilterLoads.delete(oldQuery);
                }
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
