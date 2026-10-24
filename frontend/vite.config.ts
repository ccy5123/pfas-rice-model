import { defineConfig } from "vitest/config";
import react from "@vitejs/plugin-react";

// `npm run dev` serves the UI on :5173 and proxies the API to the FastAPI backend
// (`uvicorn server.app:app --port 8000` from the repo root). `npm run build` writes
// dist/, which the backend serves itself at "/".
export default defineConfig({
  plugins: [react()],
  server: {
    proxy: { "/api": { target: "http://127.0.0.1:8000", changeOrigin: false } },
  },
  build: {
    chunkSizeWarningLimit: 5000, // plotly.js is one big, lazily-loaded chunk
  },
  test: {
    environment: "node",
  },
});
