"""Thin adapters from agent orchestration to deterministic py4siesta tools."""

from py4siesta import tool_cli


def invoke_tool(arguments):
    """Invoke a py4siesta-tool command without duplicating its implementation."""
    return tool_cli.execute(arguments)
