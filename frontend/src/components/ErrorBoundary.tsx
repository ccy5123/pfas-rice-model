import { Component, type ReactNode } from "react";

/** Keeps a rendering bug in one view from blanking the whole page (sidebar stays usable). */
export class ErrorBoundary extends Component<{ children: ReactNode; resetKey: string }, { error: Error | null }> {
  state = { error: null as Error | null };
  static getDerivedStateFromError(error: Error) { return { error }; }
  componentDidUpdate(prev: { resetKey: string }) {
    if (prev.resetKey !== this.props.resetKey && this.state.error) this.setState({ error: null });
  }
  render() {
    if (this.state.error) {
      return <div className="alert error" role="alert">Display error: {this.state.error.message}</div>;
    }
    return this.props.children;
  }
}
