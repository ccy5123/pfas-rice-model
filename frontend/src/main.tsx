import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import App from "./App";
import { loadPlotly } from "./components/Plot";
import "./styles.css";

createRoot(document.getElementById("root")!).render(
  <StrictMode>
    <App />
  </StrictMode>,
);

// start fetching the (large, lazily-split) chart library once the page is idle
const idle = (window as any).requestIdleCallback ?? ((f: () => void) => setTimeout(f, 300));
idle(() => { void loadPlotly(); });
