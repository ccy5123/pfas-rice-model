import { useEffect, useRef, useState } from "react";

/** The value, but only after it has stopped changing for `ms` (slider drags, typing). */
export function useDebounced<T>(value: T, ms = 350): T {
  const [v, setV] = useState(value);
  useEffect(() => {
    const id = setTimeout(() => setV(value), ms);
    return () => clearTimeout(id);
  }, [value, ms]);
  return v;
}

export interface Async<T> {
  data: T | null;
  /** the `key` the shown `data` was computed for (lags `key` while a new one loads) */
  dataKey: string | null;
  error: string | null;
  loading: boolean;
}

/** Run `fn` whenever `key` changes; keep showing the previous data while the new one loads
 *  (no blank flash), and drop results that arrive for a superseded key. */
export function useAsync<T>(fn: (() => Promise<T>) | null, key: string): Async<T> {
  const [state, setState] = useState<Async<T>>({ data: null, dataKey: null, error: null, loading: !!fn });
  const latest = useRef(key);
  useEffect(() => {
    latest.current = key;
    if (!fn) { setState({ data: null, dataKey: null, error: null, loading: false }); return; }
    setState((s) => ({ ...s, loading: true, error: null }));
    fn().then(
      (data) => { if (latest.current === key) setState({ data, dataKey: key, error: null, loading: false }); },
      (e: unknown) => {
        if (latest.current === key) {
          setState((s) => ({ ...s, error: e instanceof Error ? e.message : String(e), loading: false }));
        }
      },
    );
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [key]);
  return state;
}

/** Light/dark from the OS, live. */
export function useDarkMode(): boolean {
  const q = typeof window !== "undefined" ? window.matchMedia?.("(prefers-color-scheme: dark)") : null;
  const [dark, setDark] = useState(!!q?.matches);
  useEffect(() => {
    if (!q) return;
    const on = () => setDark(q.matches);
    q.addEventListener("change", on);
    return () => q.removeEventListener("change", on);
  }, [q]);
  return dark;
}
