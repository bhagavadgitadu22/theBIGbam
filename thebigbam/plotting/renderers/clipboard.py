"""Shared browser clipboard feedback for plot double-click actions."""

COPY_FEEDBACK_JS = """
function fallbackCopy(text) {
    const input = document.createElement('textarea');
    input.value = text;
    input.style.position = 'fixed';
    input.style.opacity = '0';
    document.body.appendChild(input);
    input.focus();
    input.select();
    try { document.execCommand('copy'); } catch (_) {}
    input.remove();
}
if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(name).catch(() => fallbackCopy(name));
} else {
    fallbackCopy(name);
}
const tip = document.createElement('div');
tip.textContent = 'Copied';
Object.assign(tip.style, {
    position: 'fixed', left: cb_obj.sx + 'px', top: Math.max(8, cb_obj.sy - 30) + 'px',
    transform: 'translateX(-50%)', background: '#333', color: '#fff',
    padding: '3px 8px', borderRadius: '4px', fontSize: '12px',
    pointerEvents: 'none', zIndex: '9999', opacity: '1', transition: 'opacity 0.4s ease',
});
document.body.appendChild(tip);
setTimeout(() => { tip.style.opacity = '0'; }, 700);
setTimeout(() => tip.remove(), 1100);
"""
